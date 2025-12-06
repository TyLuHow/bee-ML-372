# ApisTox Deployment Guide

Complete guide for deploying the ApisTox ML-powered bee toxicity prediction platform.

## Architecture Overview

The application uses a **split deployment** architecture:
- **Frontend (React)**: Deployed on Vercel (static build)
- **Backend (Python/FastAPI)**: Deployed on Railway (ML API server)

This separation allows:
- Fast static frontend delivery via Vercel's CDN
- Dedicated resources for ML model serving on Railway
- Independent scaling and deployment of frontend and backend

## Pre-Deployment Checklist

### Backend Files Created
- [x] `backend/api/__init__.py` - Package initialization
- [x] `backend/api/main.py` - FastAPI application (210 lines)
- [x] `backend/api/predict.py` - Prediction logic (324 lines)
- [x] `backend/api/featurizer.py` - Molecular descriptor computation (93 lines)
- [x] `backend/models/best_model_random_forest.pkl` - Trained model (3.4 MB)
- [x] `backend/models/best_model_xgboost.pkl` - Alternative model (187 KB)
- [x] `backend/models/preprocessor.pkl` - Feature scaler (2 KB)
- [x] `backend/requirements.txt` - Python dependencies
- [x] `backend/README.md` - Backend documentation

### Frontend Updates
- [x] `services/geminiService.ts` - Updated to call backend API
- [x] `types.ts` - Added `smiles` field to ChemicalData interface
- [x] `.env.local` - Added VITE_API_URL configuration
- [x] `.env.example` - Environment variable template

### Deployment Configuration
- [x] `railway.toml` - Railway deployment config (Python backend)
- [x] `nixpacks.toml` - Python runtime configuration
- [x] `README.md` - Updated with full setup instructions

## Local Testing

### 1. Test Backend Locally

```bash
# Navigate to backend directory
cd backend

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Start server
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

**Verify backend is running:**
```bash
# Health check
curl http://localhost:8000/health

# Expected response:
# {"status":"healthy","model_loaded":true,"model_name":"Random Forest"}

# Test prediction
curl http://localhost:8000/test

# Sample prediction
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "smiles": "C1=CN=C(N1)NC(=O)NCCl",
    "category": "Insecticide",
    "mw": 255.66,
    "logP": 0.57,
    "exposure": "Contact, Oral"
  }'
```

### 2. Test Frontend Locally

```bash
# From project root (with backend running on port 8000)
npm install
npm run dev
```

**Verify frontend:**
- Navigate to http://localhost:5173
- Open browser console
- Submit a compound analysis
- Check console for API calls to `http://localhost:8000/predict`
- Verify prediction results display correctly

## Railway Deployment (Backend)

### Step 1: Push to GitHub

```bash
git add .
git commit -m "Add production-ready Python backend with ML models"
git push origin main
```

**Important:** Ensure model files are included in the repository:
- `backend/models/best_model_random_forest.pkl` (3.4 MB)
- `backend/models/best_model_xgboost.pkl` (187 KB)
- `backend/models/preprocessor.pkl` (2 KB)

If files are too large for GitHub:
1. Use Git LFS: `git lfs track "*.pkl"`
2. Or upload to Railway volume storage
3. Or use cloud storage (S3, GCS) and download on startup

### Step 2: Connect Railway Project

1. Go to https://railway.app
2. Create new project
3. Select "Deploy from GitHub repo"
4. Choose your `apistox-pro` repository
5. Railway will automatically detect the Python project

### Step 3: Verify Railway Configuration

Railway will use the configuration from `railway.toml` and `nixpacks.toml`:

**railway.toml:**
```toml
[build]
builder = "nixpacks"

[deploy]
startCommand = "cd backend && uvicorn api.main:app --host 0.0.0.0 --port $PORT"
restartPolicyType = "on_failure"
restartPolicyMaxRetries = 10
```

**nixpacks.toml:**
```toml
[phases.setup]
nixPkgs = ["python311", "gcc"]

[phases.install]
cmds = ["pip install -r backend/requirements.txt"]

[start]
cmd = "cd backend && uvicorn api.main:app --host 0.0.0.0 --port $PORT"
```

### Step 4: Monitor Deployment

1. Check Railway logs for build progress
2. Look for "Loaded Random Forest model" in startup logs
3. Verify deployment is running
4. Get the Railway URL (e.g., `https://apistox-backend.railway.app`)

### Step 5: Test Railway Deployment

```bash
# Replace with your Railway URL
RAILWAY_URL="https://apistox-backend.railway.app"

# Health check
curl $RAILWAY_URL/health

# Test prediction
curl -X POST $RAILWAY_URL/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Glyphosate",
    "category": "Herbicide",
    "mw": 169.07,
    "logP": -3.4,
    "exposure": "Contact"
  }'
```

## Vercel Deployment (Frontend)

### Step 1: Update Environment Variable

In your Vercel project settings:
1. Go to Settings → Environment Variables
2. Add: `VITE_API_URL` = `https://your-railway-app.railway.app`
3. Apply to all environments (Production, Preview, Development)

### Step 2: Deploy to Vercel

Vercel will automatically deploy when you push to GitHub, or manually:

```bash
# Install Vercel CLI
npm i -g vercel

# Deploy
vercel --prod
```

### Step 3: Verify CORS

The backend CORS is configured to allow Vercel domains:
```python
allow_origins=[
    "http://localhost:3000",
    "http://localhost:5173",
    "http://localhost:4173",
    "https://*.vercel.app",
    "https://*.railway.app",
]
```

If you have issues, update `backend/api/main.py` with your specific Vercel domain.

## Testing Production Deployment

### 1. Frontend-Backend Integration

1. Open your Vercel deployment URL
2. Open browser DevTools (F12) → Network tab
3. Submit a compound analysis
4. Verify:
   - Request to Railway backend URL
   - Status 200 response
   - Prediction results display
   - No CORS errors

### 2. Test Cases

**Test 1: Known Toxic Compound**
```json
{
  "name": "Imidacloprid",
  "category": "Insecticide",
  "mw": 255.66,
  "logP": 0.57,
  "exposure": "Contact, Oral"
}
```
Expected: `toxicity: "Toxic"`, confidence > 75%

**Test 2: Known Safe Compound**
```json
{
  "name": "Glyphosate",
  "category": "Herbicide",
  "mw": 169.07,
  "logP": -3.4,
  "exposure": "Contact"
}
```
Expected: `toxicity: "Safe"`, confidence > 75%

**Test 3: With SMILES**
```json
{
  "name": "Atrazine",
  "smiles": "CCNc1nc(NC(C)C)nc(Cl)n1",
  "category": "Herbicide",
  "mw": 215.68,
  "logP": 2.61,
  "exposure": "Oral"
}
```
Expected: Valid prediction with computed molecular descriptors

**Test 4: Invalid SMILES**
```json
{
  "name": "Invalid",
  "smiles": "INVALID_SMILES_STRING",
  "category": "Insecticide",
  "mw": 200,
  "logP": 2.0,
  "exposure": "Contact"
}
```
Expected: `toxicity: "Uncertain"`, error message about invalid SMILES

### 3. Fallback Testing

Test that frontend falls back to mock predictions when backend is unavailable:

1. Stop Railway deployment temporarily
2. Try submitting analysis from frontend
3. Verify fallback message in prediction explanation
4. Restart Railway deployment

## Monitoring and Debugging

### Backend Logs (Railway)

```bash
# View Railway logs
railway logs

# Look for:
# - "Loaded Random Forest model from..."
# - "Processing prediction request for: {compound}"
# - "Prediction complete: {toxicity} (confidence: {X}%)"
```

### Frontend Console (Browser)

Check for:
- API URL being called: `${API_URL}/predict`
- Request payload structure
- Response data
- Any CORS or network errors

### Common Issues

**Issue: CORS Error**
- **Solution**: Verify backend CORS includes your Vercel domain
- Update `backend/api/main.py` `allow_origins` list

**Issue: Model not loading**
- **Solution**: Check Railway logs for import errors
- Verify all dependencies installed
- Ensure model files are in `backend/models/`

**Issue: 422 Validation Error**
- **Solution**: Check request payload matches API schema
- Verify all required fields present
- Check field names (mw vs molecular_weight)

**Issue: Slow predictions**
- **Solution**: Models are loaded once at startup (cached)
- First request may be slow (cold start)
- Subsequent requests should be fast (<500ms)

## Performance Optimization

### Backend
- Models cached in memory after first load
- Use Railway Pro for more resources if needed
- Consider model quantization for smaller file size
- Add Redis caching for repeated predictions

### Frontend
- Static build served from Vercel CDN
- API calls only when user submits analysis
- Fallback to mock predictions if API unavailable

## Security Considerations

1. **API Rate Limiting**: Consider adding rate limiting to `/predict` endpoint
2. **Input Validation**: SMILES validation prevents injection attacks
3. **CORS**: Restrict to specific frontend domains in production
4. **Error Messages**: Don't expose internal errors to frontend
5. **Model Security**: Models are read-only, loaded once at startup

## Cost Estimates

### Railway (Backend)
- **Hobby Plan**: $5/month, 500 hours
- **Starter Plan**: $10/month, unlimited hours
- Model serving ~100MB RAM, low CPU usage

### Vercel (Frontend)
- **Hobby Plan**: Free, 100GB bandwidth/month
- Static site, minimal resource usage

## Rollback Plan

If deployment issues occur:

1. **Backend Issues**: Revert Railway deployment to previous version
2. **Frontend Issues**: Revert Vercel deployment or rollback git commit
3. **Emergency**: Frontend has fallback to mock predictions

## Next Steps

1. Set up monitoring (Sentry, LogRocket)
2. Add analytics (Plausible, Google Analytics)
3. Implement caching for predictions
4. Add rate limiting
5. Set up automated tests (pytest for backend, Jest for frontend)
6. Create API documentation with Swagger UI (already included at `/docs`)

## Support

For issues:
1. Check Railway logs
2. Check browser console
3. Review backend/README.md
4. Test with `/test` endpoint
5. Verify environment variables

## Summary

The backend is production-ready with:
- 628 lines of Python code
- 3 trained ML models (3.6 MB total)
- FastAPI with CORS
- Health checks and error handling
- Comprehensive API documentation

The deployment uses:
- Railway for Python backend
- Vercel for React frontend
- Environment-based configuration
- Fallback mechanisms for reliability
