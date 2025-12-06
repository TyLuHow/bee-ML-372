<div align="center">
<img width="1200" height="475" alt="GHBanner" src="https://github.com/user-attachments/assets/0aa67016-6eaf-458a-adb2-6e31a0763ed6" />
</div>

# ApisTox - Bee Toxicity Prediction Platform

ML-powered honey bee pesticide toxicity prediction using Random Forest and XGBoost models trained on 1,035 compounds from ECOTOX, PPDB, and BPDB databases.

**View in AI Studio:** https://ai.studio/apps/drive/1jeFYUcyfxee0j8ksSN9klEpUUo2oCTA6

## Architecture

This application consists of two parts:

1. **Frontend (React + TypeScript)**: Interactive UI for compound analysis
2. **Backend (Python + FastAPI)**: ML prediction API with trained models

```
apistox-pro/
├── backend/              # Python FastAPI backend
│   ├── api/             # API endpoints and prediction logic
│   ├── models/          # Trained ML models (Random Forest, XGBoost)
│   └── requirements.txt
├── services/            # Frontend API client
├── App.tsx              # Main React application
└── types.ts             # TypeScript interfaces
```

## Run Locally

### Prerequisites

- **Node.js** (v18+) - for frontend
- **Python** (3.9+) - for backend
- **pip** - Python package manager

### Backend Setup

1. Navigate to backend directory:
   ```bash
   cd backend
   ```

2. Create virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Start the backend server:
   ```bash
   uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
   ```

   The API will be available at:
   - API: http://localhost:8000
   - Interactive docs: http://localhost:8000/docs
   - Health check: http://localhost:8000/health

### Frontend Setup

1. Navigate to project root:
   ```bash
   cd ..
   ```

2. Install dependencies:
   ```bash
   npm install
   ```

3. Configure environment variables:
   - Copy `.env.local` and set `VITE_API_URL=http://localhost:8000`
   - For production, set this to your Railway backend URL

4. Run the frontend:
   ```bash
   npm run dev
   ```

   The app will be available at http://localhost:5173

## API Documentation

### POST /predict

Predict bee toxicity for a pesticide compound.

**Request:**
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

**Response:**
```json
{
  "toxicity": "Toxic",
  "confidence": 89.3,
  "explanation": "Analysis indicates high lipophilicity...",
  "recommendation": "Recommend avoiding application during bloom periods..."
}
```

See [backend/README.md](backend/README.md) for full API documentation.

## Model Information

- **Models**: Random Forest (primary), XGBoost (alternative)
- **Training Data**: 1,035 compounds from ECOTOX, PPDB, BPDB
- **Features**: 15 molecular descriptors computed from SMILES using RDKit
- **Accuracy**: ~85-90% on test set
- **Preprocessing**: StandardScaler + SMOTE for class balance

### Molecular Descriptors

The models compute these descriptors from SMILES strings:
- MolecularWeight, LogP, TPSA
- H-bond donors/acceptors
- Rotatable bonds, aromatic rings
- Heteroatoms, ring counts
- Complexity metrics (BertzCT)
- And more...

## Deployment

### Frontend (Vercel)

The React frontend is deployed on Vercel:
1. Connect your GitHub repository
2. Set environment variable: `VITE_API_URL=<your-railway-backend-url>`
3. Deploy automatically on push

### Backend (Railway)

The Python backend is configured for Railway deployment:
1. Push to GitHub
2. Connect Railway to your repository
3. Railway will automatically detect Python and use `railway.toml` config
4. Models are included in the repository and loaded on startup

**Note:** The backend uses the configuration in `railway.toml` and `nixpacks.toml` to:
- Install Python 3.11 and dependencies
- Run FastAPI with uvicorn
- Serve on Railway's assigned PORT

## Features

- Real-time toxicity prediction using trained ML models
- SMILES-based molecular descriptor computation
- Confidence scores with scientific explanations
- Actionable safety recommendations
- Scenario-based analysis tools
- Educational resources about bee toxicity

## Development

### Frontend Development
```bash
npm run dev        # Development server
npm run build      # Production build
npm run preview    # Preview production build
```

### Backend Development
```bash
cd backend
uvicorn api.main:app --reload  # Development server with auto-reload
```

### Testing

Test the backend API:
```bash
curl http://localhost:8000/health
curl http://localhost:8000/test
```

Test a prediction:
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

### Backend not connecting
- Verify backend is running: `curl http://localhost:8000/health`
- Check VITE_API_URL in `.env.local`
- Check CORS settings in `backend/api/main.py`

### Model loading errors
- Ensure model files are in `backend/models/`
- Check Python version (3.9+ required)
- Verify all dependencies installed

### RDKit issues
- Install rdkit-pypi: `pip install rdkit-pypi`
- On some systems: `conda install -c conda-forge rdkit`

## License

Part of the ApisTox project for IME 372 Course Project.
