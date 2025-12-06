# ApisTox Backend API

Production-ready Python backend for bee toxicity prediction using trained Random Forest and XGBoost models.

## Overview

This FastAPI application provides ML-powered predictions for honey bee pesticide toxicity based on molecular features computed from compound properties.

## Architecture

```
backend/
├── api/
│   ├── __init__.py
│   ├── main.py          # FastAPI application with endpoints
│   ├── predict.py       # Prediction logic and model loading
│   └── featurizer.py    # Molecular descriptor computation (RDKit)
├── models/
│   ├── best_model_random_forest.pkl  # Trained Random Forest model
│   ├── best_model_xgboost.pkl        # Trained XGBoost model
│   └── preprocessor.pkl              # StandardScaler for feature preprocessing
└── requirements.txt     # Python dependencies
```

## Setup

### Prerequisites

- Python 3.9+
- pip or conda

### Installation

1. Create a virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install dependencies:

```bash
cd backend
pip install -r requirements.txt
```

3. Verify models are in place:

```bash
ls models/
# Should show: best_model_random_forest.pkl, best_model_xgboost.pkl, preprocessor.pkl
```

## Running the Server

### Development Mode

```bash
# From backend directory
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

Or from project root:

```bash
cd backend
python -m uvicorn api.main:app --reload
```

The API will be available at:
- API: http://localhost:8000
- Interactive docs: http://localhost:8000/docs
- Alternative docs: http://localhost:8000/redoc

### Production Mode

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8000 --workers 4
```

## API Endpoints

### GET /

Root endpoint with API information.

**Response:**
```json
{
  "name": "ApisTox Prediction API",
  "version": "1.0.0",
  "description": "ML-powered honey bee pesticide toxicity prediction",
  "endpoints": {
    "predict": "/predict",
    "health": "/health",
    "docs": "/docs"
  }
}
```

### GET /health

Health check endpoint to verify model loading.

**Response:**
```json
{
  "status": "healthy",
  "model_loaded": true,
  "model_name": "Random Forest"
}
```

### POST /predict

Main prediction endpoint.

**Request Body:**
```json
{
  "name": "Imidacloprid",
  "smiles": "C1=CN=C(N1)NC(=O)NCCl",
  "category": "Insecticide",
  "mw": 255.66,
  "logP": 0.57,
  "exposure": "Contact, Oral"
}
```

**Required Fields:**
- `name`: Compound name (string)
- `category`: Pesticide category - "Insecticide", "Fungicide", "Herbicide", or "Other"
- `exposure`: Exposure route - "Contact", "Oral", "Systemic", or combination

**Optional but Recommended:**
- `smiles`: SMILES string for accurate molecular descriptor computation
- `mw`: Molecular weight (g/mol)
- `logP`: LogP value (lipophilicity)

**Response:**
```json
{
  "toxicity": "Toxic",
  "confidence": 89.3,
  "explanation": "Imidacloprid exhibits high lipophilicity...",
  "recommendation": "Recommend avoiding application during bloom periods..."
}
```

**Response Fields:**
- `toxicity`: "Toxic", "Safe", or "Uncertain"
- `confidence`: Confidence score (0-100)
- `explanation`: Detailed scientific explanation
- `recommendation`: Actionable safety recommendations

### GET /test

Test endpoint with a sample prediction for Imidacloprid.

## Feature Engineering

The API computes 15 molecular descriptors from SMILES strings using RDKit:

1. MolecularWeight
2. LogP (lipophilicity)
3. NumHDonors (hydrogen bond donors)
4. NumHAcceptors (hydrogen bond acceptors)
5. NumRotatableBonds
6. NumAromaticRings
7. TPSA (topological polar surface area)
8. NumHeteroatoms
9. NumRings
10. NumSaturatedRings
11. NumAliphaticRings
12. FractionCSP3
13. MolarRefractivity
14. BertzCT (complexity)
15. HeavyAtomCount

These descriptors are combined with categorical features (category, exposure route) and preprocessed using the same pipeline as training.

## Model Information

- **Primary Model**: Random Forest Classifier
- **Alternative**: XGBoost Classifier
- **Training Data**: 1,035 compounds from ECOTOX, PPDB, and BPDB databases
- **Features**: 15 molecular descriptors + categorical encodings
- **Preprocessing**: StandardScaler + SMOTE for class balance
- **Performance**: ~85-90% accuracy on test set

## Error Handling

The API handles various error conditions:

- Invalid SMILES: Returns "Uncertain" with validation error message
- Missing required fields: Returns 422 validation error
- Model prediction failure: Returns "Uncertain" with error explanation
- Server errors: Returns 500 with detailed error message

## CORS Configuration

CORS is enabled for:
- localhost:3000, 5173, 4173 (local development)
- *.vercel.app (Vercel deployments)
- *.railway.app (Railway deployments)

## Deployment

### Railway

The backend is configured for Railway deployment using the root-level configuration files:

- `railway.toml`: Specifies build and start commands
- `nixpacks.toml`: Configures Python runtime

Railway will automatically:
1. Detect Python and install dependencies from `backend/requirements.txt`
2. Run the FastAPI server on the assigned PORT

### Docker (Optional)

Create a `Dockerfile` in the backend directory:

```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

Build and run:
```bash
docker build -t apistox-backend .
docker run -p 8000:8000 apistox-backend
```

## Testing

### Manual Testing

1. Start the server
2. Visit http://localhost:8000/docs
3. Try the `/test` endpoint for a quick check
4. Use the `/predict` endpoint with sample compounds

### cURL Examples

Health check:
```bash
curl http://localhost:8000/health
```

Prediction:
```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Glyphosate",
    "category": "Herbicide",
    "mw": 169.07,
    "logP": -3.4,
    "exposure": "Contact"
  }'
```

## Troubleshooting

### Models not loading
- Verify model files are in `backend/models/` directory
- Check file permissions
- Review startup logs for error messages

### RDKit import errors
- Ensure `rdkit-pypi` is installed: `pip install rdkit-pypi`
- On some systems, you may need `conda install -c conda-forge rdkit`

### Port already in use
- Change port: `uvicorn api.main:app --port 8001`
- Kill existing process: `lsof -ti:8000 | xargs kill`

## License

Part of the ApisTox project for IME 372 Course Project.
