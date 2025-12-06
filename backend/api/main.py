"""
ApisTox Backend API
FastAPI application for bee toxicity prediction
"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Optional
import logging

from .predict import get_predictor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="ApisTox Prediction API",
    description="Machine Learning API for predicting honey bee pesticide toxicity",
    version="1.0.0"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:5173",
        "http://localhost:4173",
        "https://*.vercel.app",
        "https://*.railway.app",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Request/Response Models
class CompoundInput(BaseModel):
    """Input model for toxicity prediction request"""
    name: str = Field(..., description="Compound name")
    smiles: Optional[str] = Field(None, description="SMILES string representation")
    category: str = Field(..., description="Pesticide category (Insecticide, Fungicide, Herbicide, Other)")
    molecular_weight: Optional[float] = Field(None, alias="mw", description="Molecular weight (g/mol)")
    logp: Optional[float] = Field(None, alias="logP", description="LogP (lipophilicity)")
    exposure_route: str = Field(..., alias="exposure", description="Exposure route (Contact, Oral, Systemic)")

    class Config:
        populate_by_name = True


class PredictionOutput(BaseModel):
    """Output model for toxicity prediction response"""
    toxicity: str = Field(..., description="Toxicity classification: Toxic, Safe, or Uncertain")
    confidence: float = Field(..., description="Confidence score (0-100)")
    explanation: str = Field(..., description="Detailed explanation of the prediction")
    recommendation: str = Field(..., description="Actionable recommendations")


class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    model_loaded: bool
    model_name: Optional[str] = None


# Load model on startup
@app.on_event("startup")
async def startup_event():
    """Initialize the predictor and load models"""
    try:
        predictor = get_predictor()
        logger.info(f"API started successfully. Model: {predictor.model_name}")
    except Exception as e:
        logger.error(f"Failed to load models on startup: {e}")
        raise


# API Endpoints
@app.get("/", tags=["Root"])
async def root():
    """Root endpoint - API information"""
    return {
        "name": "ApisTox Prediction API",
        "version": "1.0.0",
        "description": "ML-powered honey bee pesticide toxicity prediction",
        "endpoints": {
            "predict": "/predict",
            "health": "/health",
            "docs": "/docs"
        }
    }


@app.get("/health", response_model=HealthResponse, tags=["Health"])
async def health_check():
    """Health check endpoint"""
    try:
        predictor = get_predictor()
        return {
            "status": "healthy",
            "model_loaded": predictor.model is not None,
            "model_name": predictor.model_name
        }
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return {
            "status": "unhealthy",
            "model_loaded": False,
            "model_name": None
        }


@app.post("/predict", response_model=PredictionOutput, tags=["Prediction"])
async def predict_toxicity(compound: CompoundInput):
    """
    Predict bee toxicity for a pesticide compound.

    This endpoint accepts compound properties and returns a toxicity prediction
    with confidence score, explanation, and recommendations.

    **Required fields:**
    - name: Compound name
    - category: Pesticide category (Insecticide, Fungicide, Herbicide, Other)
    - exposure: Exposure route (Contact, Oral, Systemic)

    **Optional but recommended:**
    - smiles: SMILES string for accurate molecular descriptor computation
    - mw: Molecular weight (g/mol)
    - logP: LogP value (lipophilicity)

    **Returns:**
    - toxicity: Classification (Toxic/Safe/Uncertain)
    - confidence: Confidence score (0-100)
    - explanation: Detailed scientific explanation
    - recommendation: Actionable safety recommendations
    """
    try:
        # Get predictor instance
        predictor = get_predictor()

        # Prepare compound data
        compound_data = {
            'name': compound.name,
            'smiles': compound.smiles or '',
            'category': compound.category,
            'molecular_weight': compound.molecular_weight,
            'logp': compound.logp,
            'exposure_route': compound.exposure_route
        }

        # Generate prediction
        logger.info(f"Processing prediction request for: {compound.name}")
        result = predictor.predict(compound_data)

        logger.info(f"Prediction complete: {result['toxicity']} (confidence: {result['confidence']}%)")
        return result

    except Exception as e:
        logger.error(f"Prediction endpoint error: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"Prediction failed: {str(e)}"
        )


@app.get("/test", tags=["Testing"])
async def test_endpoint():
    """Test endpoint with a sample prediction"""
    sample_compound = {
        'name': 'Imidacloprid',
        'smiles': 'C1=CN=C(N1)NC(=O)NCCl',
        'category': 'Insecticide',
        'molecular_weight': 255.66,
        'logp': 0.57,
        'exposure_route': 'Contact, Oral'
    }

    try:
        predictor = get_predictor()
        result = predictor.predict(sample_compound)
        return {
            "test": "success",
            "sample_compound": sample_compound,
            "prediction": result
        }
    except Exception as e:
        logger.error(f"Test endpoint error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Error handlers
@app.exception_handler(404)
async def not_found_handler(request, exc):
    """Custom 404 handler"""
    return {
        "error": "Not Found",
        "message": "The requested endpoint does not exist",
        "available_endpoints": ["/", "/health", "/predict", "/docs"]
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
