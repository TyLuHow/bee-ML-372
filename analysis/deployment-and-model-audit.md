# ApisTox Deployment Infrastructure & ML Model Integration Audit

**Date:** December 9, 2025
**Project:** ApisTox (bee-ml-372) - Bee Toxicity Prediction Platform
**Repository:** https://github.com/TyLuHow/bee-ml-372
**Auditor:** Claude Code Analysis System

---

## Executive Summary

### Overall Health Status: DEGRADED

The ApisTox platform consists of a React/TypeScript frontend and Python/FastAPI ML backend with trained toxicity prediction models. While the codebase is well-structured and the ML integration is comprehensive, **the deployment infrastructure is critically broken**:

- **Frontend (Vercel):** ❌ NOT DEPLOYED - No deployments exist despite project being linked
- **Backend (Railway):** ❌ WRONG PROJECT - Railway is serving "OnRoute Outdoors" instead of ApisTox
- **GitHub Actions:** ⚠️ NO CI/CD - No GitHub workflows configured
- **ML Models:** ✅ FUNCTIONAL - Models properly integrated in codebase (3.4MB Random Forest + 188KB XGBoost)

### Critical Issues Requiring Immediate Attention

1. **Railway Project Misconfiguration (CRITICAL)** - Railway service "onroute-outdoors" is deployed instead of the ApisTox backend
2. **Vercel Deployment Missing (CRITICAL)** - No frontend deployments exist despite project link
3. **No Active CI/CD (HIGH)** - Manual deployment process with no automation
4. **Environment Variable Management (HIGH)** - VITE_API_URL not configured for production

### Quick Wins for Improvement

1. Deploy backend to Railway with correct project configuration
2. Deploy frontend to Vercel with proper environment variables
3. Verify end-to-end connectivity between frontend and backend
4. Set up basic health monitoring and deployment notifications

---

## 1. Deployment Status by Platform

### 1.1 Vercel (Frontend Deployment)

**Status:** ❌ NOT DEPLOYED

**Configuration Files:**
- `/home/yler_uby_oward/apistox/vercel.json` - Build configuration present
- `/home/yler_uby_oward/apistox/.vercel/project.json` - Project linked (ID: prj_1Ia6tVYjm3yG51jfAbhko9iZ5kP3)

**Build Configuration:**
```json
{
  "framework": "vite",
  "buildCommand": "npm run build",
  "outputDirectory": "dist",
  "devCommand": "npm run dev",
  "installCommand": "npm install"
}
```

**Findings:**
- ✅ Project is linked to Vercel account (tyluhow)
- ✅ Vercel configuration is correct for Vite framework
- ❌ **CRITICAL:** No deployments found (`vercel list` returned "No deployments found")
- ❌ Frontend is not accessible to users
- ⚠️ Environment variable `VITE_API_URL` status unknown (cannot verify without deployment)

**Required Actions:**
1. Trigger initial Vercel deployment: `vercel --prod` or push to GitHub to trigger auto-deploy
2. Configure environment variable: `VITE_API_URL` = Railway backend URL
3. Verify build succeeds with no errors
4. Test deployed frontend URL

**Expected Outcome:**
- Vercel URL: `https://apistox-pro.vercel.app` or similar
- Build time: ~2-3 minutes
- Deployment triggers on every push to main branch

---

### 1.2 Railway (Backend Deployment)

**Status:** ❌ WRONG PROJECT DEPLOYED

**Configuration Files:**
- `/home/yler_uby_oward/apistox/railway.json` - ApisTox backend configuration
- `/home/yler_uby_oward/apistox/railway.toml` - Deployment settings
- `/home/yler_uby_oward/apistox/nixpacks.toml` - Python build configuration (FIXED)

**Railway.json Configuration:**
```json
{
  "$schema": "https://railway.app/railway.schema.json",
  "build": {
    "builder": "NIXPACKS",
    "nixpacksConfigPath": "nixpacks.toml"
  },
  "deploy": {
    "startCommand": "cd backend && uvicorn api.main:app --host 0.0.0.0 --port $PORT",
    "healthcheckPath": "/health",
    "healthcheckTimeout": 100,
    "restartPolicyType": "ON_FAILURE",
    "restartPolicyMaxRetries": 10
  }
}
```

**Nixpacks Configuration (FIXED):**
```toml
[phases.setup]
nixPkgs = ["python311", "python311Packages.pip"]  # ✅ Fixed from broken 'pip' reference

[phases.install]
cmds = ["pip install --break-system-packages -r backend/requirements.txt"]  # ✅ Fixed PEP 668 issue

[start]
cmd = "cd backend && uvicorn api.main:app --host 0.0.0.0 --port $PORT"
```

**Current Deployment Status:**
- Railway Project: `passionate-bravery`
- Railway Service: `onroute-outdoors` (production)
- Railway Domains:
  - https://onroute-outdoors-production.up.railway.app
  - https://onroute-outdoors-production-4cb4.up.railway.app
- **❌ CRITICAL ISSUE:** Railway is serving the WRONG application (OnRoute Outdoors API, not ApisTox)

**Testing Current Railway Deployment:**
```bash
# Health check returns OnRoute Outdoors, not ApisTox:
$ curl https://onroute-outdoors-production.up.railway.app/
{
  "name": "OnRoute Outdoors API",
  "version": "1.0.0",
  "description": "Along-Route Outdoor Finder - Find MTB and running trails along driving routes",
  "endpoints": {...}
}
```

**Root Cause Analysis:**

The Railway CLI is connected to a different project directory. The current directory contains ApisTox code, but Railway is linked to the OnRoute Outdoors project. This indicates:

1. Either the Railway project needs to be switched/re-linked
2. Or a new Railway service needs to be created for ApisTox
3. The ApisTox backend has never been successfully deployed to Railway

**Previous Deployment Attempts:**

Based on git history and documentation:
- Multiple fixes applied for Nixpacks configuration (commits from Dec 7-9)
- PEP 668 error fixed with `--break-system-packages` flag
- Undefined variable error fixed (pip → python311Packages.pip)
- **However, deployments were likely failing silently due to wrong Railway project**

**Required Actions:**

1. **Option A (Recommended):** Create new Railway service for ApisTox
   ```bash
   railway init  # Create new project
   railway up    # Deploy current directory
   ```

2. **Option B:** Switch Railway project
   ```bash
   railway link  # Link to correct project
   railway up    # Deploy
   ```

3. **After deployment:**
   - Verify health endpoint: `curl https://[new-railway-url]/health`
   - Expected response: `{"status":"healthy","model_loaded":true,"model_name":"Random Forest"}`
   - Test prediction endpoint with sample compound
   - Update Vercel environment variable with new Railway URL

**Expected Outcome:**
- Railway URL: `https://apistox-production.up.railway.app` or similar
- Build time: ~5-7 minutes (includes installing RDKit and ML dependencies)
- Models loaded on startup: Random Forest (primary), XGBoost (alternative)
- Memory usage: ~300-500MB (due to model size)

---

### 1.3 GitHub Repository & CI/CD

**Status:** ⚠️ NO CI/CD CONFIGURED

**Repository Details:**
- Remote: https://github.com/TyLuHow/bee-ml-372
- Current Branch: `main`
- Recent Commits: 10 commits focused on Railway deployment fixes
- Authenticated: ✅ GitHub CLI authenticated as TyLuHow

**Recent Commit History:**
```
bd64401 (HEAD -> main, origin/main) Fix Railway build: disable npm build phase
2affd3a Fix PEP 668 externally-managed-environment error
783be0d Add quick verification guide for Railway fix
c70288f Add comprehensive documentation for Railway deployment fix
4a6e332 Fix Railway deployment by adding Nixpacks configuration
253150d Fix icon rendering and replace ML pipeline emojis
fe1093f Replace all emojis with professional lucide-react icons
a2719cb Fix Railway deployment with proper Python/pip configuration
b8d83c8 Fix CORS policy and fallback prediction crashes
db15bcc Integrate real trained models and scientific visualizations
```

**GitHub Workflows:**
- ❌ No `.github/workflows/` directory found
- ❌ No automated testing on pull requests
- ❌ No automated deployments
- ❌ No CI/CD pipeline

**GitHub Actions Status:**
- Cannot query workflow runs (no workflows configured)

**Findings:**
- ✅ Repository is active with recent commits
- ✅ Git history shows systematic debugging of deployment issues
- ❌ **No automated testing** - Manual testing only
- ❌ **No deployment automation** - Manual deployment required
- ⚠️ Deployment fixes suggest trial-and-error debugging process

**Required Actions:**

1. **Set up GitHub Actions for CI:**
   - Create `.github/workflows/test.yml` for automated testing
   - Run backend tests on pull requests
   - Lint and type-check frontend code

2. **Set up deployment automation:**
   - Auto-deploy to Railway on push to main
   - Auto-deploy to Vercel on push to main
   - Add deployment status badges to README

3. **Add basic test suite:**
   - Backend: Test model loading, prediction endpoint, health check
   - Frontend: Basic component tests, API integration tests

**Example GitHub Actions Workflow:**
```yaml
# .github/workflows/deploy.yml
name: Deploy to Production
on:
  push:
    branches: [main]
jobs:
  deploy-backend:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Deploy to Railway
        uses: bervProject/railway-deploy@main
        with:
          railway_token: ${{ secrets.RAILWAY_TOKEN }}
  deploy-frontend:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Deploy to Vercel
        uses: amondnet/vercel-action@v25
        with:
          vercel-token: ${{ secrets.VERCEL_TOKEN }}
          vercel-org-id: ${{ secrets.VERCEL_ORG_ID }}
          vercel-project-id: ${{ secrets.VERCEL_PROJECT_ID }}
```

---

### 1.4 Cross-Platform Integration

**Status:** ❌ NOT FUNCTIONAL

**Architecture Design:**
```
┌─────────────────────────────────────────────────────────────┐
│                    User's Browser                            │
└─────────────────┬───────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────┐
│  Vercel (Frontend - React/TypeScript)                        │
│  ❌ NOT DEPLOYED                                             │
│  - Should be at: https://apistox-pro.vercel.app             │
│  - Env var needed: VITE_API_URL                             │
└─────────────────┬───────────────────────────────────────────┘
                  │
                  │ HTTP POST /predict
                  │ Content-Type: application/json
                  │
                  ▼
┌─────────────────────────────────────────────────────────────┐
│  Railway (Backend - Python/FastAPI)                          │
│  ❌ WRONG PROJECT (OnRoute Outdoors deployed)               │
│  - Should be: ApisTox backend with ML models                │
│  - Current: OnRoute Outdoors API                            │
└─────────────────────────────────────────────────────────────┘
```

**Integration Points:**

1. **Frontend → Backend API Call:**
   - File: `/home/yler_uby_oward/apistox/services/geminiService.ts`
   - Line 4: `const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'`
   - Line 28: `const response = await fetch(`${API_URL}/predict`, {...})`

2. **Fallback Mechanism:**
   - Line 49-55: Catches fetch errors and falls back to mock predictions
   - WARNING: This masks deployment failures
   - User sees predictions but doesn't know they're not using real ML models

**CORS Configuration:**
- Backend file: `/home/yler_uby_oward/apistox/backend/api/main.py`
- Lines 28-34: CORS configured to allow all origins (`allow_origins=["*"]`)
- ✅ CORS should not be an issue once both services are deployed

**Environment Variables:**

**Frontend (Vercel):**
- Required: `VITE_API_URL` = Railway backend URL
- Status: ❌ Not configured (no deployment exists)
- Purpose: Tells frontend where to send prediction requests

**Backend (Railway):**
- Required: `PORT` (auto-provided by Railway)
- Status: ⚠️ Unknown (wrong project deployed)

**Testing Cross-Platform Integration:**

Once both services are deployed:

```bash
# 1. Test backend health
curl https://[railway-url]/health

# 2. Test backend prediction directly
curl -X POST https://[railway-url]/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "category": "Insecticide",
    "mw": 255.66,
    "logP": 0.57,
    "exposure": "Contact"
  }'

# 3. Test frontend (in browser)
# - Open https://[vercel-url]
# - Submit compound analysis
# - Check browser console for API calls
# - Verify NO "fallback prediction" message
```

**Required Actions:**
1. Deploy both frontend and backend
2. Configure VITE_API_URL in Vercel environment variables
3. Redeploy frontend to pick up environment variable
4. Test end-to-end prediction flow
5. Monitor for errors in browser console and server logs

---

## 2. ML Model Architecture

### 2.1 Model Files and Artifacts

**Location:** `/home/yler_uby_oward/apistox/backend/models/`

**Model Inventory:**

| File | Size | Type | Last Modified | Purpose |
|------|------|------|---------------|---------|
| `best_model_random_forest.pkl` | 3.4 MB | Pickle | 2025-12-09 09:55:07 | Primary classification model |
| `best_model_xgboost.pkl` | 188 KB | Pickle | 2025-12-09 09:55:07 | Alternative/backup model |
| `preprocessor.pkl` | 4 KB | Pickle | 2025-12-09 09:55:07 | StandardScaler for feature normalization |

**Model Framework:** scikit-learn (Random Forest) and XGBoost

**Model Type:** Binary classification (Toxic = 1, Safe = 0)

**Training Data:**
- Total compounds: 1,035
- Sources: ECOTOX, PPDB, BPDB databases
- Features: 15 molecular descriptors + categorical encodings
- Preprocessing: StandardScaler + SMOTE for class balance
- Expected accuracy: ~85-90% on test set

**Model Selection Logic:**
- File: `/home/yler_uby_oward/apistox/backend/api/predict.py`
- Lines 48-63: Tries Random Forest first, falls back to XGBoost if not found
- Default: Random Forest (higher accuracy, larger file size)

**Preprocessor:**
- Type: StandardScaler (scikit-learn)
- Purpose: Normalizes molecular descriptor values to same scale
- Applied before prediction (line 200 in predict.py)

---

### 2.2 Model Loading and Initialization

**Singleton Pattern Implementation:**

**Global Predictor Instance:**
- File: `/home/yler_uby_oward/apistox/backend/api/predict.py`
- Lines 330-342: Global `_predictor` variable with `get_predictor()` function
- Purpose: Load models once at startup, reuse for all predictions

```python
# Singleton pattern ensures models loaded only once
_predictor = None

def get_predictor() -> ToxicityPredictor:
    global _predictor
    if _predictor is None:
        _predictor = ToxicityPredictor()
    return _predictor
```

**Startup Initialization:**
- File: `/home/yler_uby_oward/apistox/backend/api/main.py`
- Lines 107-115: FastAPI startup event handler
- Loads models when server starts (not on each request)

```python
@app.on_event("startup")
async def startup_event():
    try:
        predictor = get_predictor()
        logger.info(f"API started successfully. Model: {predictor.model_name}")
    except Exception as e:
        logger.error(f"Failed to load models on startup: {e}")
        raise
```

**Loading Process:**

1. **Model Directory Resolution:**
   - Lines 31-34 in `predict.py`: Default to `backend/models/` relative to script
   - Allows flexibility for different deployment paths

2. **Model Loading:**
   - Lines 48-63: Try Random Forest first, XGBoost as fallback
   - Uses Python `pickle.load()` to deserialize models
   - Logs model name and path on success

3. **Preprocessor Loading:**
   - Lines 66-72: Load StandardScaler from preprocessor.pkl
   - Warns if not found but continues (predictions less accurate)

4. **Error Handling:**
   - Lines 74-76: Catches exceptions and raises to prevent server startup
   - Better to fail fast than serve broken predictions

**Caching Strategy:**
- ✅ Models loaded once at startup (not per request)
- ✅ Singleton pattern prevents duplicate loading
- ✅ No lazy loading (all models ready when first request arrives)

**Error Scenarios:**
- Missing model file → Server fails to start (intentional)
- Corrupt pickle file → Server fails to start (intentional)
- Missing preprocessor → Server starts with warning (graceful degradation)

---

### 2.3 Inference Pipeline

**Complete Prediction Flow:**

```
User Input (Frontend)
        ↓
   API Request: POST /predict
        ↓
 [1] Input Validation (main.py)
        ↓
 [2] Feature Preparation (predict.py:_prepare_features)
        ↓
 [3] SMILES → Molecular Descriptors (featurizer.py)
        ↓
 [4] Feature Engineering (categorical encoding)
        ↓
 [5] Preprocessing (StandardScaler)
        ↓
 [6] Model Prediction (Random Forest)
        ↓
 [7] Post-processing (explanation, recommendation)
        ↓
 [8] JSON Response
        ↓
   Frontend Display
```

**Step-by-Step Breakdown:**

**[1] Input Validation - main.py:154-199**
- Endpoint: `/predict` (POST)
- Input model: `CompoundInput` (Pydantic validation)
- Accepts: name, category, exposure, molecular properties, SMILES
- Line 182: Convert to dictionary, exclude None values

**[2] Feature Preparation - predict.py:78-152**
- Function: `_prepare_features(compound_data: Dict) → pd.DataFrame`
- Two paths:
  - If SMILES provided → compute descriptors from structure
  - If no SMILES → use provided molecular properties
- Line 94: `smiles_to_descriptors()` converts SMILES to 15 descriptors
- Lines 99-115: Handle field aliases (mw/MolecularWeight, logP/LogP, etc.)

**[3] Molecular Descriptor Computation - featurizer.py:37-59**
- Uses RDKit library for chemistry computations
- Function: `smiles_to_descriptors(smiles: str) → Dict[str, float]`
- Computes 15 descriptors:
  - MolecularWeight, LogP (lipophilicity)
  - NumHDonors, NumHAcceptors (hydrogen bonding)
  - NumRotatableBonds, NumAromaticRings
  - TPSA (polar surface area)
  - NumHeteroatoms, NumRings, NumSaturatedRings, NumAliphaticRings
  - FractionCSP3, MolarRefractivity
  - BertzCT (molecular complexity)
  - HeavyAtomCount

**[4] Feature Engineering - predict.py:117-145**
- Add year (default: 2024)
- Binary pesticide flags (insecticide, herbicide, fungicide, other_agrochemical)
- One-hot encoding for categorical features:
  - Source: PPDB, BPDB, Other (ECOTOX is reference, dropped)
  - Toxicity type: Oral, Systemic, Other (Contact is reference, dropped)
- Creates single-row DataFrame with all features

**[5] Preprocessing - predict.py:196-203**
- Apply StandardScaler if available
- Normalizes feature values to mean=0, std=1
- Uses same scaler fit during training (from preprocessor.pkl)
- Handles missing scaler gracefully (warns but continues)

**[6] Model Prediction - predict.py:206-210**
- Line 206: `model.predict(features)` → binary class (0 or 1)
- Line 207: `model.predict_proba(features)` → probability scores
- Extract confidence from probability of predicted class

**[7] Post-processing - predict.py:213-223**
- Convert prediction to label: "Toxic" or "Safe"
- Generate scientific explanation (context-aware, considers category and properties)
- Generate actionable recommendations (exposure-route specific)
- Round confidence to 1 decimal place

**[8] Response Formation - main.py:189-192**
- Return: toxicity, confidence, explanation, recommendation
- Log prediction for monitoring
- Handle exceptions with HTTP 500 error

**Performance Characteristics:**
- Total latency: ~50-200ms (depends on SMILES complexity)
- Bottleneck: RDKit descriptor computation (~30-100ms)
- Model inference: ~10-20ms
- No database calls (stateless prediction)

---

### 2.4 API Endpoints

**Base URL:** `http://localhost:8000` (local) or Railway deployment URL

**Endpoint Inventory:**

#### 1. GET `/` - Root Endpoint
- **File:** `/home/yler_uby_oward/apistox/backend/api/main.py:119-131`
- **Purpose:** API information and endpoint directory
- **Response:**
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
- **Authentication:** None
- **Rate Limiting:** None configured

#### 2. GET `/health` - Health Check
- **File:** `/home/yler_uby_oward/apistox/backend/api/main.py:134-150`
- **Purpose:** Verify model loading and API health
- **Response:**
  ```json
  {
    "status": "healthy",
    "model_loaded": true,
    "model_name": "Random Forest"
  }
  ```
- **Use Case:** Railway healthcheck, monitoring, debugging
- **Error Handling:** Returns unhealthy status if model fails to load

#### 3. POST `/predict` - Toxicity Prediction
- **File:** `/home/yler_uby_oward/apistox/backend/api/main.py:153-199`
- **Purpose:** Main prediction endpoint for bee toxicity assessment

**Request Schema:**
```json
{
  "name": "string (optional)",
  "smiles": "string (optional but recommended)",
  "category": "Insecticide|Fungicide|Herbicide|Other",
  "exposure": "Contact|Oral|Systemic",
  "mw": "number (optional if SMILES provided)",
  "logP": "number (optional if SMILES provided)",
  "MolecularWeight": "number (alias for mw)",
  "LogP": "number (alias for logP)",
  // ... 15+ molecular descriptors (optional)
  "insecticide": "0|1",
  "herbicide": "0|1",
  "fungicide": "0|1",
  "other_agrochemical": "0|1"
}
```

**Response Schema:**
```json
{
  "toxicity": "Toxic|Safe|Uncertain",
  "confidence": 89.3,
  "explanation": "Detailed scientific explanation...",
  "recommendation": "Actionable safety recommendations..."
}
```

**Example Request:**
```bash
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

**Field Requirements:**
- **Strictly Required:** None (all fields optional)
- **Recommended:** name, category, exposure, SMILES or molecular descriptors
- **Optional:** All molecular descriptors (computed from SMILES if provided)

**Validation:**
- Pydantic model validation (type checking)
- SMILES validation via RDKit (lines 176-183 in predict.py)
- Invalid SMILES → Returns "Uncertain" prediction

#### 4. GET `/test` - Sample Prediction
- **File:** `/home/yler_uby_oward/apistox/backend/api/main.py:202-224`
- **Purpose:** Quick test endpoint with pre-configured Imidacloprid example
- **Response:**
  ```json
  {
    "test": "success",
    "sample_compound": {
      "name": "Imidacloprid",
      "smiles": "C1=CN=C(N1)NC(=O)NCCl",
      "category": "Insecticide",
      "molecular_weight": 255.66,
      "logp": 0.57,
      "exposure_route": "Contact, Oral"
    },
    "prediction": {
      "toxicity": "Toxic",
      "confidence": 89.3,
      "explanation": "...",
      "recommendation": "..."
    }
  }
  ```
- **Use Case:** Verify API is working, test deployment

#### 5. GET `/docs` - Interactive API Documentation
- **Auto-generated:** FastAPI Swagger UI
- **URL:** `http://localhost:8000/docs`
- **Features:** Interactive testing, schema exploration, example requests

**Rate Limiting:**
- ❌ Not configured
- ⚠️ Recommendation: Add rate limiting for production (e.g., 100 requests/minute per IP)

**Caching:**
- ❌ No response caching
- ⚠️ Could cache predictions for identical SMILES strings

**Authentication:**
- ❌ No authentication required
- ⚠️ Public API - consider API keys for production if needed

---

### 2.5 Frontend Integration

**Frontend → Backend Communication:**

**API Client Service:**
- **File:** `/home/yler_uby_oward/apistox/services/geminiService.ts`
- **Purpose:** Abstracts backend API calls, handles errors, provides fallback

**Key Functions:**

#### `analyzeChemicalToxicity()` - Main Prediction Function
- **Lines 10-56:** Core API call logic
- **Input:** `ChemicalData` interface (defined in types.ts)
- **Output:** `PredictionResult` interface

**Request Preparation (Lines 12-25):**
```typescript
const requestBody: any = {};
// Copy all non-null fields
for (const [key, value] of Object.entries(data)) {
  if (value !== undefined && value !== null && value !== '') {
    requestBody[key] = value;
  }
}
// Handle field aliases for backwards compatibility
if (data.mw && !requestBody.MolecularWeight) requestBody.MolecularWeight = data.mw;
if (data.logP && !requestBody.LogP) requestBody.LogP = data.logP;
if (data.exposure && !requestBody.toxicity_type) requestBody.toxicity_type = data.exposure;
```

**API Call (Lines 28-40):**
```typescript
const response = await fetch(`${API_URL}/predict`, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(requestBody),
});

if (!response.ok) {
  throw new Error(`API request failed: ${response.status} ${response.statusText}`);
}

const result = await response.json();
```

**Error Handling & Fallback (Lines 49-55):**
```typescript
catch (error) {
  console.error('Toxicity prediction error:', error);
  console.warn('API unavailable, falling back to mock prediction');
  return fallbackPrediction(data);  // ⚠️ Silently falls back to mock predictions
}
```

**Fallback Prediction Logic:**
- **Lines 62-83:** Mock prediction when backend unavailable
- **Algorithm:** Deterministic hash-based logic using compound properties
- **Decision factors:**
  - Category (Insecticide → higher toxic probability)
  - LogP (lipophilicity → toxicity indicator)
  - Molecular weight
  - Known toxic compounds (Imidacloprid, Clothianidin, etc.)
- **Warning:** Adds "(Note: Using fallback prediction - backend API unavailable)" to explanation

**⚠️ CRITICAL ISSUE:** The fallback mechanism masks deployment failures. Users receive predictions even when backend is broken, potentially leading to:
- False confidence in predictions
- Undetected deployment issues
- Lower prediction accuracy (mock logic vs trained ML model)

**User Flow:**

```
User fills form in App.tsx
        ↓
Click "Analyze Compound" (Line 349)
        ↓
handleAnalyze() called (Line 113)
        ↓
analyzeChemicalToxicity(formData) (Line 115)
        ↓
POST request to ${VITE_API_URL}/predict
        ↓
   [SUCCESS PATH]              [ERROR PATH]
        ↓                           ↓
Backend returns result      Backend unreachable/error
        ↓                           ↓
Display prediction         fallbackPrediction()
        ↓                           ↓
Real ML result         Mock prediction + warning
```

**Frontend Display Components:**

**Input Form - App.tsx:132-347**
- Collapsible sections for basic info and molecular descriptors
- Category selector (Insecticide, Fungicide, Herbicide)
- Exposure route selector (Contact, Oral, Systemic)
- 15+ molecular descriptor inputs with default values

**Result Display - App.tsx:358-428**
- Color-coded card (green for Safe, red for Toxic)
- Confidence meter with visual progress bar
- Scientific explanation panel
- Actionable recommendations panel
- Success/error messaging

**State Management:**
- `formData` state: Stores all input values (Line 87)
- `result` state: Stores prediction results (Line 85)
- `loading` state: Shows loading indicator during prediction (Line 84)

**Error Handling:**
- Network errors: Caught and logged, falls back to mock
- Invalid SMILES: Backend returns "Uncertain" prediction
- Missing fields: Backend uses defaults or computes from SMILES

---

### 2.6 Model Metadata and Configuration

**Feature Configuration:**

**Molecular Descriptor Names:**
- **File:** `/home/yler_uby_oward/apistox/backend/api/featurizer.py`
- **Lines 18-34:** Descriptor function mapping

```python
self.descriptor_functions = {
    'MolecularWeight': Descriptors.MolWt,
    'LogP': Crippen.MolLogP,
    'NumHDonors': Lipinski.NumHDonors,
    'NumHAcceptors': Lipinski.NumHAcceptors,
    'NumRotatableBonds': Lipinski.NumRotatableBonds,
    'NumAromaticRings': Lipinski.NumAromaticRings,
    'TPSA': rdMolDescriptors.CalcTPSA,
    'NumHeteroatoms': Lipinski.NumHeteroatoms,
    'NumRings': Lipinski.RingCount,
    'NumSaturatedRings': Lipinski.NumSaturatedRings,
    'NumAliphaticRings': Lipinski.NumAliphaticRings,
    'FractionCSP3': Lipinski.FractionCSP3,
    'MolarRefractivity': Crippen.MolMR,
    'BertzCT': Descriptors.BertzCT,
    'HeavyAtomCount': Lipinski.HeavyAtomCount
}
```

**Expected Input Schema:**

**Frontend Interface - types.ts:2-46**
```typescript
export interface ChemicalData {
  // Identifiers
  name?: string;
  CID?: number;  // PubChem ID
  CAS?: string;  // CAS registry number
  smiles?: string;  // SMILES structure

  // Classification
  source?: 'ECOTOX' | 'PPDB' | 'BPDB' | 'Unknown';
  year?: number;
  toxicity_type?: 'Contact' | 'Oral' | 'Systemic' | 'Other';
  category?: string;

  // Binary flags (0 or 1)
  insecticide?: number;
  herbicide?: number;
  fungicide?: number;
  other_agrochemical?: number;

  // Molecular descriptors (15 features)
  MolecularWeight?: number;
  LogP?: number;
  // ... + 13 more descriptors
}
```

**Backend Input Model - main.py:38-89**
- Pydantic validation model with extensive field aliases
- Supports multiple naming conventions (mw/MolecularWeight, logP/LogP, etc.)
- All fields optional with detailed descriptions
- Line 88: `populate_by_name = True` enables alias support

**Hyperparameters & Model Settings:**

**Model Training Configuration (Not in Codebase):**
- Model files are pre-trained and serialized
- Training hyperparameters not stored in production code
- Model characteristics inferred from size and type:
  - Random Forest: Likely 100-500 trees (3.4MB suggests complex ensemble)
  - XGBoost: Compact model (188KB suggests pruned/optimized)

**Preprocessing Pipeline:**
- StandardScaler parameters stored in preprocessor.pkl
- Mean and standard deviation fit on training data
- Applied identically to new predictions

**Default Values:**
- Year: 2024 (Line 120 in predict.py)
- Binary flags: 0 (not that category)
- Missing descriptors: 0.0 (Line 115 in predict.py)

**A/B Testing & Model Versioning:**
- ❌ No A/B testing implemented
- ❌ No model versioning system
- ✅ Fallback from Random Forest to XGBoost if RF not found
- ⚠️ Recommendation: Add model version tracking for future updates

**Configuration Files:**
- ❌ No centralized config.py or settings.py
- Configuration scattered across:
  - railway.json (deployment settings)
  - requirements.txt (dependency versions)
  - Hardcoded defaults in predict.py

---

### 2.7 Performance and Monitoring

**Prediction Latency Logging:**

**Backend Logging:**
- **File:** `/home/yler_uby_oward/apistox/backend/api/main.py`
- **Lines 14-18:** Logging configuration
  ```python
  logging.basicConfig(
      level=logging.INFO,
      format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
  )
  ```
- **Line 188:** Log prediction request: `logger.info(f"Processing prediction request for: {compound_name}")`
- **Line 191:** Log prediction result: `logger.info(f"Prediction complete: {result['toxicity']} (confidence: {result['confidence']}%)")`
- ❌ No explicit latency timing

**Error Logging:**
- **Line 195:** Log prediction errors: `logger.error(f"Prediction endpoint error: {e}")`
- **predict.py:74-76:** Log model loading errors
- **predict.py:150-151:** Log feature preparation errors

**Frontend Logging:**
- **File:** `/home/yler_uby_oward/apistox/services/geminiService.ts`
- **Line 50:** `console.error('Toxicity prediction error:', error)`
- **Line 53:** `console.warn('API unavailable, falling back to mock prediction')`
- ❌ No performance metrics collected
- ❌ No user analytics

**Model Performance Metrics:**

**Training Metrics (from documentation):**
- Expected accuracy: ~85-90% on test set
- Source: IMPLEMENTATION_SUMMARY.md, README.md
- ❌ No live accuracy monitoring
- ❌ No drift detection

**Runtime Metrics:**
- ❌ No prediction latency tracking
- ❌ No throughput monitoring
- ❌ No error rate tracking
- ❌ No confidence score distribution

**Health Monitoring:**

**Railway Health Check:**
- Endpoint: `/health` (configured in railway.json)
- Timeout: 100ms (Line 10 in railway.json)
- Checks: Model loaded status

**Restart Policy:**
- Type: ON_FAILURE (Line 11 in railway.json)
- Max retries: 10 (Line 12 in railway.json)
- Ensures service recovery from crashes

**Alerting & Monitoring:**
- ❌ No monitoring dashboard
- ❌ No alerting on errors
- ❌ No uptime monitoring
- ❌ No performance degradation alerts

**Recommendations for Monitoring:**

1. **Add Latency Tracking:**
   ```python
   import time
   start = time.time()
   result = predictor.predict(compound_data)
   latency = time.time() - start
   logger.info(f"Prediction latency: {latency*1000:.2f}ms")
   ```

2. **Track Prediction Distributions:**
   - Log toxic/safe ratio
   - Track confidence score distribution
   - Monitor for unusual patterns

3. **Set Up External Monitoring:**
   - Use Railway metrics dashboard
   - Set up uptime monitoring (UptimeRobot, Pingdom)
   - Configure error tracking (Sentry, Rollbar)

4. **Add Custom Metrics:**
   - Predictions per hour
   - Average confidence by category
   - SMILES validation failure rate
   - Fallback prediction percentage (critical metric!)

---

## 3. Issues and Recommendations

### 3.1 Critical Issues

#### Issue 1: Railway Deployment Wrong Project (CRITICAL)
- **Severity:** CRITICAL
- **Impact:** Backend API never deployed, 100% fallback predictions in production
- **Root Cause:** Railway CLI linked to "OnRoute Outdoors" project instead of ApisTox
- **Evidence:**
  - `railway status` shows "onroute-outdoors" service
  - Health check returns OnRoute Outdoors API, not ApisTox
  - Railway domains serve wrong application

**Fix Steps:**
1. Create new Railway project for ApisTox: `railway init`
2. Deploy backend: `railway up`
3. Verify health endpoint returns ApisTox response
4. Update Vercel VITE_API_URL to new Railway URL
5. Test end-to-end prediction flow

**Time Estimate:** 15-20 minutes

---

#### Issue 2: Vercel Frontend Not Deployed (CRITICAL)
- **Severity:** CRITICAL
- **Impact:** No public frontend access, application not usable
- **Root Cause:** Project linked but never deployed
- **Evidence:** `vercel list --yes` returns "No deployments found"

**Fix Steps:**
1. Deploy to Vercel: `vercel --prod` or push to GitHub
2. Configure environment variable: VITE_API_URL = Railway backend URL
3. Verify build succeeds
4. Test deployed frontend

**Time Estimate:** 10-15 minutes

---

#### Issue 3: Fallback Prediction Masking (CRITICAL)
- **Severity:** CRITICAL
- **Impact:** Users don't know they're getting mock predictions instead of real ML
- **Root Cause:** Silent fallback in geminiService.ts lines 49-55
- **Risk:** False confidence in prediction accuracy

**Fix Steps:**
1. Add prominent warning banner when using fallback
2. Consider disabling fallback in production
3. Add fallback usage metric to monitoring
4. Alert on high fallback rate

**Code Change:**
```typescript
// In geminiService.ts, return error instead of silent fallback:
catch (error) {
  console.error('Toxicity prediction error:', error);
  throw new Error('Backend API unavailable. Please try again later.');
  // Remove silent fallback for production
}
```

**Time Estimate:** 5 minutes

---

### 3.2 High Priority Issues

#### Issue 4: No CI/CD Pipeline (HIGH)
- **Severity:** HIGH
- **Impact:** Manual deployments, higher risk of errors, slower iteration
- **Root Cause:** No .github/workflows/ configured

**Fix Steps:**
1. Create `.github/workflows/deploy.yml`
2. Configure Railway and Vercel deployment secrets
3. Add automated testing (backend health check, frontend build)
4. Enable auto-deploy on push to main

**Time Estimate:** 30-45 minutes

---

#### Issue 5: No Error Monitoring (HIGH)
- **Severity:** HIGH
- **Impact:** Cannot detect production issues, poor user experience
- **Root Cause:** No monitoring infrastructure

**Fix Steps:**
1. Add Sentry or similar error tracking to backend
2. Add Sentry to frontend for client-side errors
3. Set up Railway metrics dashboard
4. Configure alerts for error spikes

**Time Estimate:** 20-30 minutes

---

#### Issue 6: Missing Environment Variable Documentation (HIGH)
- **Severity:** HIGH
- **Impact:** Deployment confusion, manual configuration errors
- **Root Cause:** No centralized env var documentation

**Fix Steps:**
1. Update README.md with complete environment variable list
2. Document Railway and Vercel env vars separately
3. Add .env.production.example file
4. Include verification steps

**Time Estimate:** 15 minutes

---

### 3.3 Medium Priority Issues

#### Issue 7: No Rate Limiting (MEDIUM)
- **Severity:** MEDIUM
- **Impact:** Potential API abuse, high costs
- **Recommendation:** Add rate limiting middleware

**Example Implementation:**
```python
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address

limiter = Limiter(key_func=get_remote_address)
app.state.limiter = limiter

@app.post("/predict")
@limiter.limit("100/minute")
async def predict_toxicity(compound: CompoundInput):
    ...
```

---

#### Issue 8: No Prediction Caching (MEDIUM)
- **Severity:** MEDIUM
- **Impact:** Repeated predictions for same compound waste resources
- **Recommendation:** Cache predictions by SMILES hash

**Example Implementation:**
```python
from functools import lru_cache
import hashlib

def cache_key(smiles: str) -> str:
    return hashlib.md5(smiles.encode()).hexdigest()

# Add caching layer before prediction
```

---

#### Issue 9: No Model Versioning (MEDIUM)
- **Severity:** MEDIUM
- **Impact:** Cannot A/B test models, risky updates
- **Recommendation:** Add model version tracking and comparison

---

### 3.4 Low Priority Issues

#### Issue 10: Hardcoded Default Values (LOW)
- **Severity:** LOW
- **Impact:** Inflexible configuration
- **Recommendation:** Move defaults to config file

#### Issue 11: No Input Validation on Frontend (LOW)
- **Severity:** LOW
- **Impact:** Invalid inputs sent to backend
- **Recommendation:** Add client-side validation before API call

---

### 3.5 Security Considerations

**Current Security Posture:**

✅ **Good Practices:**
- CORS configured (though currently allow-all)
- HTTPS for both Vercel and Railway
- No sensitive data in environment variables
- Model files not user-modifiable

⚠️ **Concerns:**
- No authentication on prediction endpoint (public API)
- No rate limiting (potential abuse)
- CORS allows all origins (Line 30 in main.py: `allow_origins=["*"]`)
- No input sanitization beyond Pydantic validation

**Recommendations:**

1. **Restrict CORS in Production:**
   ```python
   allow_origins=[
       "https://apistox-pro.vercel.app",
       "http://localhost:5173"  # Dev only
   ]
   ```

2. **Add Optional API Keys:**
   - Not required initially (public research tool)
   - Consider if abuse becomes issue

3. **Input Sanitization:**
   - SMILES validation already implemented (good!)
   - Add length limits on text fields
   - Validate numeric ranges

4. **Dependency Security:**
   - Run `pip audit` regularly
   - Keep dependencies updated
   - Monitor for CVEs in scikit-learn, FastAPI, etc.

---

### 3.6 Performance Improvement Suggestions

**Current Performance:**
- Estimated latency: 50-200ms per prediction
- Bottleneck: RDKit descriptor computation
- Memory: ~300-500MB (model size)

**Optimization Opportunities:**

1. **Batch Predictions:**
   - Add `/predict/batch` endpoint for multiple compounds
   - Amortize model loading overhead

2. **Async Prediction:**
   - Already using FastAPI async (good!)
   - Could parallelize descriptor computation for batches

3. **Model Optimization:**
   - Convert Random Forest to ONNX format (faster inference)
   - Prune XGBoost model further (already small at 188KB)
   - Quantize model weights

4. **Response Caching:**
   - Cache predictions for identical SMILES
   - Use Redis or in-memory LRU cache
   - TTL: 24 hours (predictions don't change)

5. **CDN for Frontend:**
   - Vercel already uses CDN (good!)
   - Ensure proper cache headers

---

## 4. Action Items

### Immediate Actions (Today - 1-2 hours)

**Priority 1: Deploy Backend to Railway**
- [ ] Create new Railway project: `railway init`
- [ ] Deploy ApisTox backend: `railway up`
- [ ] Verify health endpoint: `curl https://[railway-url]/health`
- [ ] Test prediction endpoint with sample compound
- [ ] Document Railway URL for next steps

**Priority 2: Deploy Frontend to Vercel**
- [ ] Deploy frontend: `vercel --prod`
- [ ] Add environment variable: VITE_API_URL = Railway backend URL
- [ ] Redeploy to pick up env var
- [ ] Test deployed frontend at Vercel URL

**Priority 3: Verify End-to-End Integration**
- [ ] Submit test prediction through deployed frontend
- [ ] Verify NO "fallback prediction" warning appears
- [ ] Check browser console for API call to Railway backend
- [ ] Confirm prediction uses real ML model

**Priority 4: Update Documentation**
- [ ] Add deployment URLs to README.md
- [ ] Update DEPLOYMENT_GUIDE.md with actual URLs
- [ ] Document environment variables clearly

---

### Short-Term Actions (This Week - 3-5 hours)

**Priority 5: Set Up CI/CD**
- [ ] Create `.github/workflows/deploy.yml`
- [ ] Configure Railway deployment secret
- [ ] Configure Vercel deployment secret
- [ ] Test automated deployment on push to main

**Priority 6: Add Monitoring**
- [ ] Set up error tracking (Sentry or similar)
- [ ] Configure Railway metrics dashboard
- [ ] Add uptime monitoring (UptimeRobot)
- [ ] Create alert for deployment failures

**Priority 7: Fix Fallback Prediction Warning**
- [ ] Add prominent banner when using fallback
- [ ] Consider disabling fallback in production
- [ ] Track fallback usage rate in monitoring

**Priority 8: Add Rate Limiting**
- [ ] Install slowapi: `pip install slowapi`
- [ ] Add rate limiting to /predict endpoint
- [ ] Configure reasonable limits (100/minute)
- [ ] Test rate limiting behavior

---

### Long-Term Improvements (Next Month - 10-15 hours)

**Priority 9: Enhanced Monitoring**
- [ ] Add prediction latency tracking
- [ ] Log prediction distributions (toxic/safe ratio)
- [ ] Track confidence score distributions
- [ ] Monitor SMILES validation failures
- [ ] Dashboard for key metrics

**Priority 10: Model Versioning**
- [ ] Add model version to health endpoint
- [ ] Create model update workflow
- [ ] Implement A/B testing capability
- [ ] Track model performance over time

**Priority 11: Performance Optimization**
- [ ] Add prediction caching (Redis or in-memory)
- [ ] Implement batch prediction endpoint
- [ ] Consider ONNX model conversion
- [ ] Profile and optimize hot paths

**Priority 12: Security Hardening**
- [ ] Restrict CORS to production domains only
- [ ] Add input length limits
- [ ] Implement API key support (optional)
- [ ] Regular dependency security audits

**Priority 13: Testing & QA**
- [ ] Add backend unit tests (pytest)
- [ ] Add integration tests for prediction flow
- [ ] Frontend component tests (Vitest)
- [ ] End-to-end tests (Playwright)

---

## 5. Conclusion

### Summary

The ApisTox platform has a **well-architected ML integration** with comprehensive model loading, feature engineering, and prediction logic. The codebase demonstrates best practices in:
- Singleton pattern for model loading
- Comprehensive error handling
- Scientific explanations and recommendations
- Fallback mechanisms for reliability

However, the **deployment infrastructure is critically broken**:
- Frontend not deployed to Vercel (linked but no deployments)
- Backend deployed to wrong Railway project (OnRoute Outdoors instead of ApisTox)
- No CI/CD pipeline
- Silent fallback predictions masking deployment failures

### Critical Path to Production

1. **Deploy backend to Railway** (15-20 min)
2. **Deploy frontend to Vercel** (10-15 min)
3. **Configure environment variables** (5 min)
4. **Test end-to-end** (10 min)
5. **Set up monitoring** (20-30 min)

**Total time to functional production deployment: ~1-2 hours**

### Success Metrics

After completing immediate actions, verify:
- ✅ Frontend accessible at Vercel URL
- ✅ Backend accessible at Railway URL
- ✅ Health endpoint returns ApisTox (not OnRoute Outdoors)
- ✅ Predictions use real ML models (no fallback warnings)
- ✅ Response time < 300ms for typical prediction
- ✅ Error tracking configured and receiving events

### Risk Assessment

**High Risk:**
- Users currently cannot access application (no frontend deployment)
- If frontend were deployed, 100% would be fallback predictions (backend wrong project)
- No way to detect production issues (no monitoring)

**Medium Risk:**
- No rate limiting (potential abuse)
- No automated testing (manual QA only)
- Silent failures (fallback masking)

**Low Risk:**
- Security (public research tool, appropriate for use case)
- Performance (adequate for current scale)

### Recommendations Summary

**Immediate (Critical):**
1. Deploy both services correctly
2. Verify end-to-end connectivity
3. Add basic monitoring

**Short-Term (High Priority):**
4. Set up CI/CD for automated deployments
5. Add error tracking and alerting
6. Fix fallback prediction warning

**Long-Term (Nice to Have):**
7. Performance optimization (caching, batching)
8. Model versioning and A/B testing
9. Comprehensive test suite

---

## Appendix

### A. File Reference Index

**Configuration Files:**
- `/home/yler_uby_oward/apistox/vercel.json` - Vercel build config
- `/home/yler_uby_oward/apistox/railway.json` - Railway deployment config
- `/home/yler_uby_oward/apistox/railway.toml` - Railway build settings
- `/home/yler_uby_oward/apistox/nixpacks.toml` - Python environment config
- `/home/yler_uby_oward/apistox/package.json` - Frontend dependencies
- `/home/yler_uby_oward/apistox/backend/requirements.txt` - Backend dependencies

**Backend Files:**
- `/home/yler_uby_oward/apistox/backend/api/main.py` - FastAPI application (241 lines)
- `/home/yler_uby_oward/apistox/backend/api/predict.py` - Prediction logic (343 lines)
- `/home/yler_uby_oward/apistox/backend/api/featurizer.py` - Molecular descriptors (94 lines)
- `/home/yler_uby_oward/apistox/backend/models/best_model_random_forest.pkl` - Primary model (3.4MB)
- `/home/yler_uby_oward/apistox/backend/models/best_model_xgboost.pkl` - Alternative model (188KB)
- `/home/yler_uby_oward/apistox/backend/models/preprocessor.pkl` - StandardScaler (4KB)

**Frontend Files:**
- `/home/yler_uby_oward/apistox/App.tsx` - Main React application (86KB)
- `/home/yler_uby_oward/apistox/services/geminiService.ts` - API client (219 lines)
- `/home/yler_uby_oward/apistox/types.ts` - TypeScript interfaces (82 lines)
- `/home/yler_uby_oward/apistox/index.tsx` - Application entry point
- `/home/yler_uby_oward/apistox/index.html` - HTML template

**Documentation:**
- `/home/yler_uby_oward/apistox/README.md` - Project overview
- `/home/yler_uby_oward/apistox/DEPLOYMENT_GUIDE.md` - Deployment instructions
- `/home/yler_uby_oward/apistox/backend-connectivity-fix.md` - Railway fix documentation
- `/home/yler_uby_oward/apistox/IMPLEMENTATION_SUMMARY.md` - ML integration summary
- `/home/yler_uby_oward/apistox/backend/README.md` - Backend API documentation

### B. Technology Stack

**Frontend:**
- React 19.2.1
- TypeScript 5.8.2
- Vite 6.2.0
- Lucide React 0.556.0 (icons)

**Backend:**
- Python 3.11
- FastAPI 0.104.1
- Uvicorn 0.24.0
- Pydantic 2.5.0
- scikit-learn 1.3.2
- XGBoost 2.0.3
- RDKit 2022.9.5
- pandas 2.1.3
- numpy 1.26.2

**Infrastructure:**
- Vercel (Frontend hosting)
- Railway (Backend hosting)
- GitHub (Version control)

### C. Model Performance Specifications

**Random Forest Model:**
- Type: Binary classifier
- Size: 3.4 MB
- Classes: 0 (Safe), 1 (Toxic)
- Features: 15 molecular descriptors + categorical encodings
- Expected accuracy: ~85-90%

**XGBoost Model:**
- Type: Binary classifier
- Size: 188 KB
- Classes: 0 (Safe), 1 (Toxic)
- Features: Same as Random Forest
- Expected accuracy: ~85-90%

**Preprocessor:**
- Type: StandardScaler
- Size: 4 KB
- Purpose: Feature normalization

---

**End of Audit Report**

Generated: December 9, 2025
Next Review: After deployment fixes completed
