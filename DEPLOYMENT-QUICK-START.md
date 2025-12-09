# Railway Deployment Quick Start Guide

## What Was Fixed

All deployment issues have been resolved:
- Added gcc/g++ for Python package compilation
- Created diagnostic startup script
- Increased healthcheck timeout to 300s
- Fixed PORT variable handling
- Committed Dockerfile to git

## Current Status

**Commit:** cd60d3f
**Pushed:** December 9, 2024 at 19:16 UTC
**Status:** Waiting for Railway to build and deploy (5-10 minutes)

## Monitor Deployment

### Option 1: Railway Dashboard (Recommended)
1. Go to https://railway.app/dashboard
2. Open project "diligent-surprise"
3. Click on the "web" service
4. Watch the "Deployments" tab for build progress
5. Check "Logs" tab for startup script output

### Option 2: Command Line
```bash
# Check if deployed (run every 60 seconds)
curl https://web-production-f014a.up.railway.app/health

# When working, you should see:
# {"status":"healthy","model_loaded":true,"model_name":"Random Forest"}
```

## What to Look For in Logs

The startup script will output:
```
================================
ApisTox Backend Starting
================================
Python version: Python 3.11.x
Working directory: /app
PORT environment variable: 8080

Checking Python packages...
  uvicorn: 0.24.0
  fastapi: 0.104.1
  scikit-learn: 1.3.2
  xgboost: 2.0.3
  rdkit: 2022.9.5

Checking model files...
total 3.6M
-rw-r--r-- 1 root root 3.4M best_model_random_forest.pkl
-rw-r--r-- 1 root root 187K best_model_xgboost.pkl
-rw-r--r-- 1 root root 2.0K preprocessor.pkl

Starting uvicorn server on 0.0.0.0:8080
================================
INFO:     Started server process [1]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:8080
```

## Test Endpoints

### 1. Health Check
```bash
curl https://web-production-f014a.up.railway.app/health
```

### 2. Root Endpoint
```bash
curl https://web-production-f014a.up.railway.app/
```

### 3. Test Prediction (Built-in Sample)
```bash
curl https://web-production-f014a.up.railway.app/test
```

### 4. Real Prediction
```bash
curl -X POST https://web-production-f014a.up.railway.app/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "smiles": "C1=CN=C(N1)NC(=O)NCCl",
    "category": "Insecticide",
    "molecular_weight": 255.66,
    "logp": 0.57,
    "exposure_route": "Contact"
  }'
```

### 5. Interactive API Docs
Open in browser: https://web-production-f014a.up.railway.app/docs

## Troubleshooting

### If Deployment Fails

1. **Check Build Logs:**
   - Look for package installation errors
   - Verify all system dependencies installed
   - Check if model files are included

2. **Check Runtime Logs:**
   - Look for ERROR messages in startup script
   - Verify all packages imported successfully
   - Check model files are accessible

3. **Common Issues:**
   - Missing packages: Look for "ERROR: package not found"
   - Model loading failure: Check for pickle errors
   - Port binding: Verify PORT variable is correct

### If Still Returning 502

Railway may still be building. Wait 10 minutes total from push time, then:

```bash
# Force redeploy if needed
railway redeploy --service web
```

## Files Changed

- `/home/yler_uby_oward/apistox/Dockerfile` (CREATED)
- `/home/yler_uby_oward/apistox/railway.json` (MODIFIED)

## Full Report

See `deployment-fix-report.md` for comprehensive details on:
- Root cause analysis
- All changes made
- Architecture overview
- Performance considerations
- Lessons learned
- Next steps

## Support

If deployment fails after 10 minutes:
1. Check Railway logs for specific error messages
2. Review the full `deployment-fix-report.md`
3. Verify git push succeeded: `git log --oneline -1` should show cd60d3f
4. Confirm Railway is connected to GitHub repository
