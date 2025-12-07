# Backend API Connectivity Fix - ApisTox

**Date:** December 7, 2025
**Issue:** Predictions showing "Note: Using fallback prediction - backend API unavailable" instead of real ML model predictions
**Status:** FIXED - Deployment configuration updated, awaiting Railway rebuild

---

## Executive Summary

The backend API was **never successfully deployed** to Railway due to a **Nixpacks build configuration error**. The frontend fallback logic was masking this complete failure, making it appear as a simple connectivity issue when in reality the backend service never started.

**Root Cause:** Invalid Nixpacks configuration trying to use undefined `pip` package instead of `python311Packages.pip`

**Impact:** 100% of predictions were fallback predictions; no ML model predictions were ever served to production users

**Fix Applied:** Updated `nixpacks.toml` with correct Python package references and pushed to trigger new Railway deployment

---

## Root Cause Analysis

### The Problem

Railway uses Nixpacks for building Python applications. The original `nixpacks.toml` configuration contained:

```toml
[phases.setup]
nixPkgs = ['python311', 'pip']  # ❌ ERROR: 'pip' is undefined
```

### Why This Failed

In the Nix package ecosystem:
- `python311` is a valid package
- `pip` is **NOT** a standalone package
- `pip` is included as `python311Packages.pip`

### The Error

Railway build logs showed:
```
error: undefined variable 'pip'
at /app/.nixpacks/nixpkgs-ffeebf0acf3ae8b29f8c7049cd911b9636efd7e7.nix:19:9
```

This caused:
1. Build to fail immediately
2. No Python environment created
3. No dependencies installed
4. No backend service started
5. No ML models loaded
6. 100% frontend fallback predictions

---

## Project Structure Discovery

During diagnosis, we discovered the project had been reorganized:

### Old Structure (app/backend/)
- `/app/backend/main.py` - FastAPI application
- `/outputs/models/` - ML model files
- `requirements-production.txt` - Production dependencies

### New Structure (backend/)
- `/backend/api/main.py` - Refactored FastAPI application
- `/backend/models/` - Relocated model files
- `/backend/requirements.txt` - Backend-specific dependencies

The frontend still uses `/app/frontend/` structure.

---

## The Fix

### 1. Updated `nixpacks.toml`

**Location:** `/nixpacks.toml` (project root)

**Changes:**
```toml
# BEFORE (broken):
[phases.setup]
nixPkgs = ['python311', 'pip']  # ❌ pip undefined

# AFTER (fixed):
[phases.setup]
nixPkgs = ["python311", "python311Packages.pip"]  # ✅ Correct Nix package path
```

**Full corrected configuration:**
```toml
# Nixpacks configuration for Railway deployment
# This file configures the Python environment and build process
# FIX: Use python311Packages.pip instead of just 'pip' to avoid undefined variable error

[phases.setup]
nixPkgs = ["python311", "python311Packages.pip"]

[phases.install]
cmds = ["pip install -r backend/requirements.txt"]

[start]
cmd = "cd backend && uvicorn api.main:app --host 0.0.0.0 --port $PORT"
```

### 2. Updated `railway.json`

**Location:** `/railway.json` (project root)

**Changes:**
- Removed redundant `buildCommand` that was conflicting with `nixpacks.toml`
- Kept proper configuration referencing nixpacks.toml

**Configuration:**
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

### 3. Git Commit and Push

**Commit:** `4a6e332`
**Message:** "Fix Railway deployment by adding Nixpacks configuration"

**Files Changed:**
- `nixpacks.toml` - Created with corrected configuration
- `railway.json` - Updated to reference nixpacks.toml properly

**Pushed to:** `main` branch on `https://github.com/TyLuHow/bee-ML-372.git`

---

## Expected Build Process

After the fix, Railway will:

1. **Setup Phase** - Install Python 3.11 and pip from Nix packages ✅
2. **Install Phase** - Run `pip install -r backend/requirements.txt` ✅
3. **Build Phase** - (No build steps needed for this project) ✅
4. **Start Phase** - Launch FastAPI with `uvicorn api.main:app` ✅

### Expected Dependencies Installed

From `/backend/requirements.txt`:
- FastAPI 0.104.1
- Uvicorn 0.24.0
- Pydantic 2.5.0
- scikit-learn 1.3.2
- xgboost 2.0.3
- rdkit-pypi 2022.9.5
- pandas 2.1.3
- numpy 1.26.2
- imbalanced-learn 0.11.0

### Expected Models Loaded

From `/backend/models/`:
- `best_model_random_forest.pkl` (3.4 MB) - Primary prediction model
- `best_model_xgboost.pkl` (187 KB) - Alternative model
- `preprocessor.pkl` (2 KB) - Feature scaler

---

## Verification Steps

### Step 1: Check Railway Deployment Status

1. Go to Railway dashboard: https://railway.app
2. Navigate to the `apis_tox_dataset` or `bee-ML-372` project
3. Check deployment logs for:
   ```
   [OK] Python environment initialized
   [OK] Installing dependencies...
   [OK] Dependencies installed successfully
   [OK] Starting FastAPI application...
   [OK] Loaded Random Forest model from backend/models/best_model_random_forest.pkl
   [OK] Application startup complete
   ```

### Step 2: Get Railway Deployment URL

In Railway dashboard:
1. Go to Settings → Domains
2. Copy the generated Railway domain (e.g., `https://apis-tox-dataset-production.up.railway.app`)

### Step 3: Test Backend Health Endpoint

Replace `<RAILWAY_URL>` with your actual Railway URL:

```bash
# Health check
curl https://<RAILWAY_URL>/health

# Expected response:
{
  "status": "healthy",
  "model_loaded": true,
  "model_name": "Random Forest",
  "timestamp": "2025-12-07T..."
}
```

### Step 4: Test Prediction Endpoint

```bash
# Test prediction with known toxic compound (Imidacloprid)
curl -X POST https://<RAILWAY_URL>/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "smiles": "C1=CN=C(N1)NC(=O)NCCl",
    "category": "Insecticide",
    "mw": 255.66,
    "logP": 0.57,
    "exposure": "Contact, Oral"
  }'

# Expected response:
{
  "toxicity": "Toxic",
  "confidence": 0.85,
  "explanation": "...",
  "model_used": "Random Forest"
}
```

### Step 5: Update Frontend Environment Variable

**Critical:** The frontend needs to know the Railway backend URL.

#### Option A: Vercel Dashboard
1. Go to Vercel dashboard
2. Select the `apistox-pro` project
3. Navigate to Settings → Environment Variables
4. Update or add:
   - **Key:** `VITE_API_URL`
   - **Value:** `https://<RAILWAY_URL>` (your actual Railway URL, no trailing slash)
   - **Environments:** Production, Preview, Development
5. Redeploy the frontend

#### Option B: Local Testing First
1. Create `/app/frontend/.env.local`:
   ```
   VITE_API_URL=https://<RAILWAY_URL>
   ```
2. Run locally: `cd app/frontend && npm run dev`
3. Test predictions in browser
4. If successful, update Vercel environment variables

### Step 6: Test End-to-End Predictions

#### Test Case 1: Toxic Insecticide
**Input:**
- Name: Imidacloprid
- Category: Insecticide
- Molecular Weight: 255.66
- LogP: 0.57
- Exposure: Contact, Oral

**Expected:**
- Toxicity: "Toxic"
- Confidence: >75%
- NO fallback message
- Scientific explanation from ML model

#### Test Case 2: Safer Herbicide
**Input:**
- Name: Glyphosate
- Category: Herbicide
- Molecular Weight: 169.07
- LogP: -3.4
- Exposure: Contact

**Expected:**
- Toxicity: "Safe" or "Low Toxicity"
- Confidence: >70%
- NO fallback message
- Scientific explanation from ML model

#### Test Case 3: With SMILES String
**Input:**
- Name: Atrazine
- SMILES: CCNc1nc(NC(C)C)nc(Cl)n1
- Category: Herbicide
- Molecular Weight: 215.68
- LogP: 2.61
- Exposure: Oral

**Expected:**
- Valid prediction using computed molecular descriptors
- Confidence score from model
- NO fallback message

### Step 7: Verify in Browser Console

1. Open deployed frontend in browser
2. Open DevTools (F12) → Network tab
3. Submit a compound analysis
4. Check:
   - ✅ Request goes to Railway URL (not localhost)
   - ✅ Status: 200 OK
   - ✅ Response contains ML model prediction
   - ✅ No CORS errors
   - ✅ No "fallback prediction" message in UI

---

## Success Criteria

| Criterion | Status | How to Verify |
|-----------|--------|---------------|
| Railway backend builds successfully | ⏳ Pending | Check Railway build logs for "Build successful" |
| Backend /health endpoint returns healthy | ⏳ Pending | `curl https://<RAILWAY_URL>/health` |
| Model files loaded correctly | ⏳ Pending | Check startup logs for "Loaded Random Forest model" |
| /predict endpoint accepts requests | ⏳ Pending | Test with curl command above |
| Frontend connects to Railway backend | ⏳ Pending | Update VITE_API_URL and test |
| Predictions use real ML models | ⏳ Pending | Verify no "fallback prediction" message |
| Toxic compounds predicted correctly | ⏳ Pending | Test Imidacloprid → "Toxic" |
| Non-toxic compounds predicted correctly | ⏳ Pending | Test Glyphosate → "Safe" |
| SMILES-based predictions work | ⏳ Pending | Test with SMILES string input |
| No CORS errors in browser | ⏳ Pending | Check browser console |

---

## Monitoring and Troubleshooting

### Railway Logs

**View logs:**
```bash
# If Railway CLI is linked to project:
railway logs

# Or check in Railway dashboard → Deployments → View Logs
```

**Look for success indicators:**
```
✅ "Python 3.11 environment initialized"
✅ "Successfully installed fastapi uvicorn scikit-learn..."
✅ "Loaded Random Forest model from backend/models/best_model_random_forest.pkl"
✅ "Application startup complete"
✅ "Uvicorn running on http://0.0.0.0:$PORT"
```

**Look for error indicators:**
```
❌ "ModuleNotFoundError" - Missing dependencies
❌ "FileNotFoundError: backend/models/*.pkl" - Model files missing
❌ "ImportError" - Package compatibility issue
❌ "Port already in use" - Port binding issue
```

### Common Issues and Solutions

#### Issue 1: Build Still Failing
**Symptoms:** Railway build logs show errors during pip install

**Solutions:**
- Check `backend/requirements.txt` for incompatible versions
- Verify all model files are committed to git
- Check if rdkit-pypi needs different version for Python 3.11
- Try specifying Python version in nixpacks.toml: `nixPkgs = ["python311"]`

#### Issue 2: Model Files Not Found
**Symptoms:** Runtime error "FileNotFoundError: backend/models/..."

**Solutions:**
- Verify model files are in git: `git ls-files | grep backend/models`
- Check if files exceed GitHub size limit (>100MB requires Git LFS)
- Ensure .gitignore doesn't exclude .pkl files
- Alternative: Upload models to Railway volume storage

#### Issue 3: CORS Errors in Browser
**Symptoms:** Browser console shows "CORS policy blocked"

**Solutions:**
- Update `backend/api/main.py` CORS middleware:
  ```python
  app.add_middleware(
      CORSMiddleware,
      allow_origins=[
          "https://apistox-pro.vercel.app",  # Add your Vercel domain
          "https://*.vercel.app",
      ],
      allow_credentials=True,
      allow_methods=["*"],
      allow_headers=["*"],
  )
  ```
- Commit and push changes
- Wait for Railway to redeploy

#### Issue 4: Frontend Still Using Fallback
**Symptoms:** Predictions show "Note: Using fallback prediction..."

**Solutions:**
1. Check `VITE_API_URL` is set in Vercel
2. Verify environment variable has no trailing slash
3. Check browser Network tab - is request going to Railway URL?
4. If request is going to localhost, environment variable not loaded
5. Redeploy frontend after updating environment variable

#### Issue 5: Slow First Prediction
**Symptoms:** First prediction takes >10 seconds

**This is normal:**
- Railway may cold-start the container
- Models are loaded into memory on first use
- Subsequent predictions should be <500ms
- Consider Railway "Always On" to prevent cold starts

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────┐
│  USER BROWSER                                           │
│  https://apistox-pro.vercel.app                         │
└────────────────────┬────────────────────────────────────┘
                     │
                     │ HTTPS Request
                     │ POST /predict
                     ↓
┌─────────────────────────────────────────────────────────┐
│  VERCEL (Frontend)                                      │
│  ├─ React + TypeScript                                  │
│  ├─ services/api.ts (VITE_API_URL)                      │
│  └─ Fallback logic if backend unavailable               │
└────────────────────┬────────────────────────────────────┘
                     │
                     │ HTTP Request
                     │ ${VITE_API_URL}/predict
                     ↓
┌─────────────────────────────────────────────────────────┐
│  RAILWAY (Backend)                                      │
│  https://<project>.up.railway.app                       │
│  ├─ Python 3.11 + FastAPI                               │
│  ├─ backend/api/main.py                                 │
│  ├─ backend/api/predict.py                              │
│  └─ backend/models/                                     │
│      ├─ best_model_random_forest.pkl (3.4 MB)           │
│      ├─ best_model_xgboost.pkl (187 KB)                 │
│      └─ preprocessor.pkl (2 KB)                         │
└─────────────────────────────────────────────────────────┘
```

---

## Next Steps

### Immediate Actions Required

1. **Monitor Railway Deployment** (Est: 5-10 minutes)
   - Check Railway dashboard for build completion
   - Verify no build errors in logs
   - Confirm service is running

2. **Get Railway URL** (Est: 1 minute)
   - Copy deployment URL from Railway dashboard
   - Note URL format: `https://<project>-production.up.railway.app`

3. **Test Backend Directly** (Est: 2 minutes)
   - Run health check curl command
   - Run prediction curl command
   - Verify responses contain model predictions

4. **Update Frontend Environment Variable** (Est: 3 minutes)
   - Add `VITE_API_URL` to Vercel environment variables
   - Set to Railway URL (no trailing slash)
   - Redeploy frontend

5. **End-to-End Testing** (Est: 10 minutes)
   - Test all 3 test cases listed above
   - Verify no fallback messages
   - Check browser console for errors
   - Confirm predictions match expected model behavior

### Optional Enhancements

6. **Set Up Monitoring** (Est: 30 minutes)
   - Add Sentry for error tracking
   - Configure Railway metrics and alerts
   - Set up uptime monitoring (UptimeRobot, Pingdom)

7. **Performance Optimization** (Est: 1 hour)
   - Enable Railway "Always On" to prevent cold starts
   - Add Redis caching for repeated predictions
   - Implement request rate limiting

8. **Documentation Updates** (Est: 30 minutes)
   - Update README.md with Railway URL
   - Document environment variables
   - Add troubleshooting section for common issues

---

## Files Modified

### Created
- `/nixpacks.toml` - Nixpacks build configuration (corrected)

### Modified
- `/railway.json` - Railway deployment configuration (cleaned up)

### Next to Update (After Deployment Success)
- `/app/frontend/.env.local` - Add Railway URL for local testing
- `/.env.example` - Update with Railway URL example
- `/README.md` - Add deployment status and URLs
- `/DEPLOYMENT_GUIDE.md` - Add actual Railway URL

---

## Lessons Learned

### What Went Wrong

1. **Silent Failure:** Frontend fallback logic masked complete backend failure
2. **Nixpacks Syntax:** Used incorrect Nix package reference (`pip` vs `python311Packages.pip`)
3. **No Monitoring:** No alerts when backend failed to deploy
4. **Incomplete Testing:** Backend was never actually tested in production

### Best Practices Going Forward

1. **Health Checks:** Frontend should show backend status (connected/disconnected)
2. **Monitoring:** Set up alerts for deployment failures and health check failures
3. **Documentation:** Keep Railway URL documented and version-controlled
4. **Testing:** Test backend independently before frontend integration
5. **Error Visibility:** Show clear error messages when backend is unavailable

---

## Summary

**Problem:** Backend never successfully deployed due to Nixpacks configuration error using undefined `pip` variable

**Solution:** Updated `nixpacks.toml` to use correct Nix package path `python311Packages.pip`

**Status:** Fix committed and pushed (commit 4a6e332), awaiting Railway rebuild

**Next:** Monitor Railway deployment, test health endpoint, update frontend environment variable

**Timeline:**
- Issue existed: Since initial Railway deployment
- Issue diagnosed: December 7, 2025
- Fix implemented: December 7, 2025
- Fix deployed: Pending Railway rebuild (~5-10 minutes)
- End-to-end testing: Pending

---

## Contact and Support

**Repository:** https://github.com/TyLuHow/bee-ML-372.git
**Railway Project:** apis_tox_dataset (or bee-ML-372)
**Vercel Project:** apistox-pro

**For Issues:**
1. Check Railway deployment logs first
2. Test backend health endpoint
3. Verify frontend environment variables
4. Review browser console for errors
5. Consult this document's troubleshooting section

---

**Document Version:** 1.0
**Last Updated:** December 7, 2025
**Author:** Claude Code (Sonnet 4.5)
