<objective>
Build a production-ready Python backend API that integrates the trained bee toxicity prediction model and connects it to the existing React frontend. Replace the mock service with actual ML-powered predictions using the pre-trained Random Forest and XGBoost models.
</objective>

<context>
The ApisTox application currently uses a mock prediction service in the React frontend. The actual machine learning models have already been trained and saved:
- Random Forest model: `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/outputs/models/best_model_random_forest.pkl`
- XGBoost model: `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/outputs/models/best_model_xgboost.pkl`
- Preprocessor: `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/outputs/preprocessors/preprocessor.pkl`

Training pipeline reference: `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/notebooks/ApisTox_Full_Pipeline.ipynb`
Dataset: `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/outputs/dataset_final.csv`

The models predict bee toxicity based on molecular features computed from compound properties (SMILES, molecular weight, LogP, category, etc.). The frontend needs to send compound data to an API endpoint and receive toxicity predictions with confidence scores and explanations.

Current deployment: React app on Vercel and Railway (static frontend).
</context>

<requirements>
1. **Backend API (Python)**:
   - Create a FastAPI or Flask backend service
   - Load the trained models and preprocessor on startup
   - Implement `/predict` endpoint that accepts compound properties
   - Compute molecular descriptors from SMILES using RDKit (same as training pipeline)
   - Run predictions using the loaded model
   - Return toxicity classification, confidence score, and explanation
   - Handle errors gracefully (invalid SMILES, missing features, etc.)
   - Enable CORS for React frontend communication

2. **Model Integration**:
   - Copy trained model files to backend directory
   - Load models and preprocessor at startup (not per request)
   - Use the same feature engineering pipeline as training (MolecularFeaturizer from notebook)
   - Apply preprocessing (scaling, encoding) before prediction
   - Generate confidence scores and interpretable explanations

3. **Frontend Updates**:
   - Update `./services/geminiService.ts` to call the backend API instead of mock data
   - Add environment variable for API URL (configurable for local dev vs production)
   - Maintain the same UI/UX and response structure
   - Handle API errors and loading states

4. **Deployment Configuration**:
   - Update Railway deployment to run the Python backend (not static serving)
   - Create requirements.txt with all dependencies
   - Update railway.toml and nixpacks.toml for Python runtime
   - Keep Vercel for React frontend (static build)
   - Configure environment variables for API URL on both platforms

5. **Documentation**:
   - Update README with backend setup instructions
   - Document API endpoint specifications
   - Include local development instructions
</requirements>

<implementation>
**Backend Structure:**
```
./backend/
├── api/
│   ├── __init__.py
│   ├── main.py              # FastAPI app
│   ├── predict.py           # Prediction endpoint logic
│   └── featurizer.py        # MolecularFeaturizer (from notebook)
├── models/
│   ├── best_model_random_forest.pkl
│   ├── best_model_xgboost.pkl
│   └── preprocessor.pkl
├── requirements.txt
├── Dockerfile (optional)
└── README.md
```

**Implementation Steps:**
1. Create backend directory structure
2. Copy trained models from `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/outputs/` to `./backend/models/`
3. Implement MolecularFeaturizer class (port from notebook cell 7)
4. Create FastAPI application with `/predict` endpoint
5. Load models on startup with caching
6. Implement prediction logic with proper preprocessing
7. Add CORS middleware for frontend communication
8. Create requirements.txt with dependencies: fastapi, uvicorn, scikit-learn, xgboost, rdkit-pypi, pandas, numpy
9. Update frontend service to call API endpoint
10. Configure Railway for Python backend deployment
11. Test end-to-end integration

**Prediction Flow:**
1. Frontend sends: `{name, smiles, category, molecular_weight, logp, exposure_route}`
2. Backend computes molecular descriptors from SMILES
3. Backend preprocesses features (same as training)
4. Backend runs model prediction
5. Backend returns: `{toxicity, confidence, explanation, recommendation}`

**Error Handling:**
- Invalid SMILES → return error with message
- Missing required fields → validation error
- Model prediction failure → fallback to uncertainty response
- Detailed logging for debugging
</implementation>

<output>
Create the following files:
- `./backend/api/main.py` - FastAPI application with CORS
- `./backend/api/predict.py` - Prediction endpoint logic
- `./backend/api/featurizer.py` - Molecular descriptor computation
- `./backend/requirements.txt` - Python dependencies
- `./backend/README.md` - Backend setup and API documentation
- `./backend/models/` - Directory with copied model files
- `./services/geminiService.ts` - Updated to call backend API
- `./.env.local` - Add API_URL environment variable
- `./railway.toml` - Updated for Python backend
- `./nixpacks.toml` - Updated for Python runtime
- `./README.md` - Updated with backend setup instructions
</output>

<verification>
Before declaring complete:
1. Test backend locally: `uvicorn backend.api.main:app --reload`
2. Verify `/predict` endpoint accepts compound data and returns predictions
3. Test with multiple compounds (toxic and non-toxic examples)
4. Confirm frontend connects to backend API
5. Verify CORS is properly configured
6. Test error handling (invalid SMILES, missing fields)
7. Check Railway deployment configuration is correct
8. Ensure model files are properly loaded (not too large for git)
9. Verify predictions match expected behavior from training notebook
</verification>

<success_criteria>
- Backend API running and accepting requests
- `/predict` endpoint returns realistic toxicity predictions
- Frontend successfully calls backend and displays results
- Model files properly integrated and loaded
- CORS configured for frontend-backend communication
- Railway configured for Python backend deployment
- requirements.txt includes all necessary dependencies
- Error handling for invalid inputs
- Documentation complete with API specifications
- End-to-end prediction flow working
</success_criteria>
