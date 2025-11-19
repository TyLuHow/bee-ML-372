#!/usr/bin/env python3
"""
FastAPI Backend for Honey Bee Toxicity Prediction
==================================================

REST API for serving ML model predictions with:
- Prediction endpoint with toxicity classification
- SHAP explanation endpoint
- Model information endpoint
- Prediction history tracking
- Data explorer endpoints
- CORS support for frontend

Author: IME 372 Project Team
Date: November 2025
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List, Dict
import pandas as pd
from datetime import datetime
import os

# Import explorer router and utilities
from app.backend import explorer
from app.backend.utils import load_model, load_json, save_json, format_error, validate_input_features, truncate_history

# Initialize FastAPI app
app = FastAPI(
    title="Honey Bee Toxicity Prediction API",
    description="ML-powered API for predicting pesticide toxicity to honey bees",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify exact origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include explorer router
app.include_router(explorer.router, prefix="/api/explorer", tags=["explorer"])

# Load model and preprocessor at startup
MODEL_PATH = "outputs/models/best_model_xgboost.pkl"
PREPROCESSOR_PATH = "outputs/preprocessors/preprocessor.pkl"
RESULTS_PATH = "outputs/metrics/training_results.json"
HISTORY_FILE = "app/backend/prediction_history.json"

model = None
preprocessor = None
model_info = {}
prediction_history = []

@app.on_event("startup")
async def load_models():
    """Load model and preprocessor on startup."""
    global model, preprocessor, model_info, prediction_history

    try:
        # Load model using utility function
        model = load_model(MODEL_PATH, "XGBoost Model")

        # Fix pickle module path issue (preprocessing -> src.preprocessing)
        import sys
        import src.preprocessing as preprocessing_module
        sys.modules['preprocessing'] = preprocessing_module

        # Load preprocessor
        preprocessor = load_model(PREPROCESSOR_PATH, "Preprocessor")

        # Verify preprocessor has required attributes
        if not hasattr(preprocessor, 'scaler'):
            print("⚠️ Warning: Preprocessor missing 'scaler' attribute")
            # Try to handle dict format for backward compatibility
            if isinstance(preprocessor, dict):
                print("   Detected dict format, attempting conversion...")
                from src.preprocessing import DataPreprocessor
                new_preprocessor = DataPreprocessor(random_state=preprocessor.get('random_state', 42))
                new_preprocessor.scaler = preprocessor.get('scaler')
                new_preprocessor.label_encoders = preprocessor.get('label_encoders', {})
                new_preprocessor.feature_selector = preprocessor.get('feature_selector')
                new_preprocessor.selected_features = preprocessor.get('selected_features')
                preprocessor = new_preprocessor
                print("   ✓ Converted dict to DataPreprocessor instance")

        # Load model info
        if os.path.exists(RESULTS_PATH):
            model_info = load_json(RESULTS_PATH, "Model metrics").get('xgboost', {})

        # Load prediction history
        if os.path.exists(HISTORY_FILE):
            prediction_history = load_json(HISTORY_FILE, "Prediction history")

        print("✓ API Ready!")

    except Exception as e:
        print(f"Error loading model: {e}")
        import traceback
        traceback.print_exc()
        raise


# Pydantic models for request/response validation
class PredictionInputSMILES(BaseModel):
    """Input for SMILES-based prediction."""
    smiles: str = Field(..., description="SMILES string of molecule", min_length=1)
    year: int = Field(2024, description="Publication year", ge=1800, le=2030)
    herbicide: int = Field(0, description="Is herbicide (0 or 1)", ge=0, le=1)
    fungicide: int = Field(0, description="Is fungicide (0 or 1)", ge=0, le=1)
    insecticide: int = Field(0, description="Is insecticide (0 or 1)", ge=0, le=1)
    other_agrochemical: int = Field(0, description="Is other agrochemical (0 or 1)", ge=0, le=1)
    source: str = Field("PPDB", description="Data source (ECOTOX or PPDB)")
    toxicity_type: str = Field("Contact", description="Toxicity type (Contact or Oral)")
    
    class Config:
        schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "year": 2024,
                "herbicide": 0,
                "fungicide": 0,
                "insecticide": 1,
                "other_agrochemical": 0,
                "source": "PPDB",
                "toxicity_type": "Contact"
            }
        }


class PredictionInput(BaseModel):
    """Input features for prediction."""
    source: str = Field(..., description="Data source (ECOTOX, PPDB, BPDB)")
    year: int = Field(..., description="Publication year", ge=1800, le=2030)
    toxicity_type: str = Field(..., description="Toxicity type (Contact, Oral, Other)")
    herbicide: int = Field(..., description="Is herbicide (0 or 1)", ge=0, le=1)
    fungicide: int = Field(..., description="Is fungicide (0 or 1)", ge=0, le=1)
    insecticide: int = Field(..., description="Is insecticide (0 or 1)", ge=0, le=1)
    other_agrochemical: int = Field(..., description="Is other agrochemical (0 or 1)", ge=0, le=1)
    MolecularWeight: float = Field(..., description="Molecular weight", ge=0)
    LogP: float = Field(..., description="Partition coefficient (lipophilicity)")
    NumHDonors: int = Field(..., description="Number of hydrogen bond donors", ge=0)
    NumHAcceptors: int = Field(..., description="Number of hydrogen bond acceptors", ge=0)
    NumRotatableBonds: int = Field(..., description="Number of rotatable bonds", ge=0)
    NumAromaticRings: int = Field(..., description="Number of aromatic rings", ge=0)
    TPSA: float = Field(..., description="Topological polar surface area", ge=0)
    NumHeteroatoms: int = Field(..., description="Number of heteroatoms", ge=0)
    NumRings: int = Field(..., description="Number of rings", ge=0)
    NumSaturatedRings: int = Field(..., description="Number of saturated rings", ge=0)
    NumAliphaticRings: int = Field(..., description="Number of aliphatic rings", ge=0)
    FractionCSP3: float = Field(..., description="Fraction of sp3 carbons", ge=0, le=1)
    MolarRefractivity: float = Field(..., description="Molar refractivity", ge=0)
    BertzCT: float = Field(..., description="Bertz molecular complexity", ge=0)
    HeavyAtomCount: int = Field(..., description="Number of heavy atoms", ge=0)
    
    class Config:
        schema_extra = {
            "example": {
                "source": "PPDB",
                "year": 2020,
                "toxicity_type": "Contact",
                "herbicide": 0,
                "fungicide": 0,
                "insecticide": 1,
                "other_agrochemical": 0,
                "MolecularWeight": 350.0,
                "LogP": 3.5,
                "NumHDonors": 2,
                "NumHAcceptors": 4,
                "NumRotatableBonds": 5,
                "NumAromaticRings": 1,
                "TPSA": 70.0,
                "NumHeteroatoms": 5,
                "NumRings": 2,
                "NumSaturatedRings": 0,
                "NumAliphaticRings": 1,
                "FractionCSP3": 0.4,
                "MolarRefractivity": 95.0,
                "BertzCT": 500.0,
                "HeavyAtomCount": 25
            }
        }


class PredictionOutput(BaseModel):
    """Prediction output with confidence scores."""
    prediction: int = Field(..., description="Predicted class (0=non-toxic, 1=toxic)")
    prediction_label: str = Field(..., description="Human-readable prediction")
    confidence: float = Field(..., description="Prediction confidence (0-1)")
    probabilities: Dict[str, float] = Field(..., description="Class probabilities")
    timestamp: str = Field(..., description="Prediction timestamp")
    
    
class ModelInfo(BaseModel):
    """Model metadata and performance."""
    model_type: str
    features: List[str]
    performance: Dict[str, float]
    training_date: str


# API Endpoints

@app.get("/")
async def root():
    """Root endpoint with API information."""
    return {
        "message": "Honey Bee Toxicity Prediction API",
        "version": "1.0.0",
        "endpoints": {
            "/predict": "Make toxicity predictions",
            "/predict/smiles": "Predict from SMILES string",
            "/model/info": "Get model information",
            "/feature/importance": "Get feature importance",
            "/history": "View prediction history",
            "/analysis/toxicophores": "Get toxicophore analysis",
            "/recommend/alternatives/{cid}": "Get safer alternatives",
            "/health": "Health check"
        },
        "explorer_endpoints": {
            "/api/explorer/overview": "Dataset overview statistics",
            "/api/explorer/molecular-diversity": "Molecular descriptor distributions",
            "/api/explorer/toxicity-by-class": "Toxicity rates by chemical type",
            "/api/explorer/temporal-trends": "Toxicity trends over time",
            "/api/explorer/chemical-space": "PCA and t-SNE coordinates",
            "/api/explorer/toxicophores": "Toxicophore enrichment data",
            "/api/explorer/correlations": "Feature correlation matrix",
            "/api/explorer/property-distributions": "2D property relationships"
        }
    }


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "model_loaded": model is not None,
        "preprocessor_loaded": preprocessor is not None,
        "timestamp": datetime.now().isoformat()
    }


@app.post("/predict", response_model=PredictionOutput)
async def predict(input_data: PredictionInput):
    """
    Make a toxicity prediction for a pesticide compound.
    
    Args:
        input_data: Pesticide features
        
    Returns:
        Prediction with confidence scores
    """
    if model is None or preprocessor is None:
        raise HTTPException(status_code=503, detail="Model not loaded")
    
    try:
        # Convert input to dataframe
        input_dict = input_data.dict()
        input_df = pd.DataFrame([input_dict])

        # Encode categorical features (matching training preprocessing)
        input_df = pd.get_dummies(input_df, columns=['source', 'toxicity_type'], drop_first=True)

        # Expected feature columns
        expected_cols = [
            'year', 'herbicide', 'fungicide', 'insecticide', 'other_agrochemical',
            'MolecularWeight', 'LogP', 'NumHDonors', 'NumHAcceptors',
            'NumRotatableBonds', 'NumAromaticRings', 'TPSA', 'NumHeteroatoms',
            'NumRings', 'NumSaturatedRings', 'NumAliphaticRings', 'FractionCSP3',
            'MolarRefractivity', 'BertzCT', 'HeavyAtomCount',
            'source_ECOTOX', 'source_PPDB', 'toxicity_type_Oral', 'toxicity_type_Other'
        ]

        # Validate and align features using utility function
        input_df = validate_input_features(input_df, expected_cols, add_missing=True)
        
        # Scale features
        input_scaled = preprocessor.scaler.transform(input_df)
        
        # Make prediction
        prediction = model.predict(input_scaled)[0]
        probabilities = model.predict_proba(input_scaled)[0]
        
        # Create response
        response = {
            "prediction": int(prediction),
            "prediction_label": "Toxic" if prediction == 1 else "Non-toxic",
            "confidence": float(probabilities[prediction]),
            "probabilities": {
                "non_toxic": float(probabilities[0]),
                "toxic": float(probabilities[1])
            },
            "timestamp": datetime.now().isoformat()
        }
        
        # Save to history
        history_entry = {**input_dict, **response}
        prediction_history.append(history_entry)

        # Truncate history using utility function
        prediction_history[:] = truncate_history(prediction_history, max_length=100)

        # Save history to file
        save_json(prediction_history, HISTORY_FILE, "Prediction history")
        
        return response
        
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Prediction error: {str(e)}")


@app.post("/predict/smiles", response_model=PredictionOutput)
async def predict_from_smiles(input_data: PredictionInputSMILES):
    """
    Predict toxicity from SMILES string.
    
    This endpoint converts a SMILES string to molecular descriptors using RDKit,
    then makes a prediction using the trained model.
    
    Args:
        input_data: SMILES string and pesticide metadata
        
    Returns:
        Prediction with confidence scores
    """
    if model is None or preprocessor is None:
        raise HTTPException(status_code=503, detail="Model not loaded")
    
    try:
        # Import featurizer
        from src.molecular_features import MolecularFeaturizer
        
        # Initialize featurizer
        featurizer = MolecularFeaturizer()
        
        # Convert SMILES to descriptors
        descriptors = featurizer.smiles_to_descriptors(input_data.smiles)
        
        if descriptors is None:
            raise HTTPException(
                status_code=400, 
                detail=f"Invalid SMILES string: {input_data.smiles}. Please provide a valid SMILES."
            )
        
        # Create full input by combining descriptors with metadata
        full_input = PredictionInput(
            **descriptors,
            year=input_data.year,
            herbicide=input_data.herbicide,
            fungicide=input_data.fungicide,
            insecticide=input_data.insecticide,
            other_agrochemical=input_data.other_agrochemical,
            source=input_data.source,
            toxicity_type=input_data.toxicity_type
        )
        
        # Use existing prediction logic
        result = await predict(full_input)
        
        # Add SMILES to response
        result_dict = result.dict() if hasattr(result, 'dict') else result
        result_dict['smiles'] = input_data.smiles
        result_dict['descriptors_calculated'] = list(descriptors.keys())
        
        return result_dict
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"SMILES prediction error: {str(e)}")


@app.get("/model/info", response_model=ModelInfo)
async def get_model_info():
    """Get model information and performance metrics."""
    if model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")
    
    feature_names = [
        'year', 'herbicide', 'fungicide', 'insecticide', 'other_agrochemical',
        'MolecularWeight', 'LogP', 'NumHDonors', 'NumHAcceptors', 
        'NumRotatableBonds', 'NumAromaticRings', 'TPSA', 'NumHeteroatoms',
        'NumRings', 'NumSaturatedRings', 'NumAliphaticRings', 'FractionCSP3',
        'MolarRefractivity', 'BertzCT', 'HeavyAtomCount',
        'source_ECOTOX', 'source_PPDB', 'toxicity_type_Oral', 'toxicity_type_Other'
    ]
    
    performance = model_info.get('val_metrics', {
        'accuracy': 0.8558,
        'f1': 0.7368,
        'roc_auc': 0.8788
    })
    
    # Filter out non-numeric metrics (confusion_matrix, classification_report)
    performance_numeric = {k: v for k, v in performance.items() 
                          if isinstance(v, (int, float))}
    
    return {
        "model_type": "XGBoost Classifier",
        "features": feature_names,
        "performance": performance_numeric,
        "training_date": model_info.get('timestamp', datetime.now().isoformat())
    }


@app.get("/history")
async def get_history(limit: int = 10):
    """Get recent prediction history."""
    return {
        "total_predictions": len(prediction_history),
        "recent_predictions": prediction_history[-limit:]
    }


@app.get("/feature/importance")
async def get_feature_importance():
    """Get global feature importance from SHAP analysis."""
    importance_path = "outputs/metrics/feature_importance_shap.csv"
    
    if not os.path.exists(importance_path):
        raise HTTPException(status_code=404, detail="Feature importance data not found")
    
    try:
        importance_df = pd.read_csv(importance_path)
        importance_data = importance_df.head(15).to_dict(orient='records')
        
        return {
            "top_features": importance_data,
            "description": "Features ranked by mean absolute SHAP value"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading importance data: {str(e)}")


@app.get("/analysis/toxicophores")
async def get_toxicophore_analysis():
    """
    Get pre-computed toxicophore analysis results.
    
    Returns statistical analysis of structural alerts associated with toxicity.
    """
    results_path = "outputs/analysis/toxicophore_results.json"
    
    if not os.path.exists(results_path):
        raise HTTPException(
            status_code=404, 
            detail="Toxicophore analysis not found. Run: python src/toxicophores.py"
        )
    
    try:
        with open(results_path, 'r') as f:
            results = json.load(f)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading toxicophore data: {str(e)}")


@app.post("/analysis/toxicophores/molecule")
async def analyze_molecule_toxicophores(smiles: str):
    """
    Identify toxicophores in a specific molecule.
    
    Args:
        smiles: SMILES string of molecule
        
    Returns:
        List of toxicophores found in the molecule
    """
    try:
        from src.toxicophores import ToxicophoreAnalyzer
        
        analyzer = ToxicophoreAnalyzer()
        toxicophores = analyzer.find_toxicophores(smiles)
        
        if not any(toxicophores.values()):
            return {
                "smiles": smiles,
                "toxicophores_found": [],
                "count": 0
            }
        
        found = [
            {
                "name": name,
                "display_name": analyzer.toxicophore_names[name]
            }
            for name, present in toxicophores.items() if present
        ]
        
        return {
            "smiles": smiles,
            "toxicophores_found": found,
            "count": len(found)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error analyzing toxicophores: {str(e)}")


@app.get("/recommend/alternatives/{compound_cid}")
async def get_alternatives(compound_cid: str, n_alternatives: int = 5):
    """
    Find safer alternative compounds similar to a given compound.
    
    Args:
        compound_cid: CID of compound to find alternatives for
        n_alternatives: Number of alternatives to return (default: 5)
        
    Returns:
        List of similar non-toxic compounds
    """
    alternatives_path = "outputs/analysis/alternatives.csv"
    
    if not os.path.exists(alternatives_path):
        raise HTTPException(
            status_code=404,
            detail="Alternatives database not found. Run: python src/recommendations.py"
        )
    
    try:
        alternatives_df = pd.read_csv(alternatives_path)
        
        # Filter by original compound ID
        compound_alts = alternatives_df[alternatives_df['original_id'] == compound_cid]
        
        if compound_alts.empty:
            raise HTTPException(
                status_code=404,
                detail=f"No alternatives found for compound {compound_cid}"
            )
        
        # Sort by rank and limit
        compound_alts = compound_alts.sort_values('rank').head(n_alternatives)
        
        result = {
            "original_id": compound_cid,
            "original_name": compound_alts.iloc[0]['original_name'],
            "n_alternatives": len(compound_alts),
            "alternatives": []
        }
        
        for idx, row in compound_alts.iterrows():
            result["alternatives"].append({
                "rank": int(row['rank']),
                "name": row['name'],
                "cid": row['CID'],
                "similarity_score": float(row['similarity_score']),
                "similarity_distance": float(row['similarity_distance'])
            })
        
        return result
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving alternatives: {str(e)}")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)

