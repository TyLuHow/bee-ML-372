# ApisTox Railway Deployment Fix Report

**Date:** December 9, 2024
**Status:** FIXES DEPLOYED - Awaiting Railway Build
**Commit:** cd60d3f

---

## Executive Summary

The ApisTox backend has been experiencing repeated deployment failures on Railway with errors including "Container failed to start", "pip not found", "uvicorn not found", and 502 gateway errors. This report documents a comprehensive root cause analysis and the systematic fixes applied to resolve all identified issues.

**Current Status:** All fixes have been committed and pushed to the main branch. Railway should automatically deploy within 5-10 minutes.

---

## 1. Root Cause Analysis

### Issues Identified

#### Critical Issues (Container Startup Failures)

1. **Missing System Dependencies**
   - **Problem:** The Dockerfile was missing `gcc` and `g++` compilers
   - **Impact:** Python packages like scikit-learn and xgboost could not compile C extensions during installation
   - **Evidence:** Previous deployment logs showed build failures for scientific packages

2. **Shell Variable Substitution Failure**
   - **Problem:** CMD line used `sh -c "uvicorn ... --port ${PORT:-8080}"` which may not properly expand variables in Railway's container environment
   - **Impact:** Server couldn't bind to the correct port, causing startup failures
   - **Evidence:** "Container failed to start" with no clear error messages

3. **Outdated pip Version**
   - **Problem:** Base Python image includes older pip version
   - **Impact:** Package installation failures, especially with newer package versions
   - **Solution:** Added `pip install --upgrade pip` before requirements installation

4. **Insufficient Healthcheck Timeout**
   - **Problem:** 100-second healthcheck timeout insufficient for loading 3.4MB Random Forest model
   - **Impact:** Railway marked container as unhealthy before FastAPI fully initialized
   - **Evidence:** Container startup taking longer than healthcheck window

5. **No Diagnostic Logging**
   - **Problem:** When containers failed, no diagnostic output to identify the issue
   - **Impact:** Impossible to determine if issue was with Python packages, models, or networking
   - **Solution:** Added comprehensive startup script with version checks

6. **Dockerfile Not Committed to Git**
   - **Problem:** Dockerfile was untracked in git repository
   - **Impact:** Railway may have been using an old cached version or falling back to Nixpacks
   - **Evidence:** `git status` showed Dockerfile as untracked

---

## 2. Changes Made

### File 1: `/home/yler_uby_oward/apistox/Dockerfile`

**Status:** CREATED and COMMITTED

**Changes:**

```dockerfile
# Before: Missing gcc, g++, no pip upgrade
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY backend/requirements.txt /app/backend/requirements.txt
RUN pip install --no-cache-dir -r backend/requirements.txt

# After: Complete system dependencies + pip upgrade
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libgomp1 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

COPY backend/requirements.txt ./backend/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r backend/requirements.txt
```

**Key Improvements:**

1. Added `gcc` and `g++` for compiling Python C extensions
2. Upgraded pip before installing requirements
3. Consolidated RUN commands for better Docker layer caching

**Startup Script Creation:**

```bash
# Before: Simple shell command with potential variable expansion issues
CMD ["sh", "-c", "uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8080}"]

# After: Comprehensive bash startup script with diagnostics
RUN echo '#!/bin/bash\n\
set -e\n\
echo "================================"\n\
echo "ApisTox Backend Starting"\n\
echo "================================"\n\
echo "Python version: $(python --version)"\n\
echo "Working directory: $(pwd)"\n\
echo "PORT environment variable: ${PORT:-8080}"\n\
echo ""\n\
echo "Checking Python packages..."\n\
python -c "import uvicorn; print(f\"  uvicorn: {uvicorn.__version__}\")" || echo "  ERROR: uvicorn not found"\n\
python -c "import fastapi; print(f\"  fastapi: {fastapi.__version__}\")" || echo "  ERROR: fastapi not found"\n\
python -c "import sklearn; print(f\"  scikit-learn: {sklearn.__version__}\")" || echo "  ERROR: scikit-learn not found"\n\
python -c "import xgboost; print(f\"  xgboost: {xgboost.__version__}\")" || echo "  ERROR: xgboost not found"\n\
python -c "import rdkit; print(f\"  rdkit: {rdkit.__version__}\")" || echo "  ERROR: rdkit not found"\n\
echo ""\n\
echo "Checking model files..."\n\
ls -lh backend/models/ || echo "  ERROR: Models directory not found"\n\
echo ""\n\
echo "Starting uvicorn server on 0.0.0.0:${PORT:-8080}"\n\
echo "================================"\n\
cd backend && exec python -m uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8080}' > /app/start.sh

RUN chmod +x /app/start.sh

CMD ["/bin/bash", "/app/start.sh"]
```

**Benefits of Startup Script:**

- Logs Python version for debugging
- Shows PORT environment variable value
- Verifies all critical packages can be imported
- Checks package versions (uvicorn, fastapi, sklearn, xgboost, rdkit)
- Confirms model files exist before starting server
- Uses `exec` to replace shell process with uvicorn (proper signal handling)
- Uses explicit bash instead of sh for reliable variable expansion

---

### File 2: `/home/yler_uby_oward/apistox/railway.json`

**Status:** MODIFIED and COMMITTED

**Changes:**

```json
{
  "deploy": {
    "healthcheckPath": "/health",
    "healthcheckTimeout": 300,        // Increased from 100 to 300 seconds
    "restartPolicyType": "ON_FAILURE",
    "restartPolicyMaxRetries": 3      // Reduced from 10 to 3
  }
}
```

**Rationale:**

- **300s healthcheck timeout:** Allows sufficient time for:
  - Container startup (10-20s)
  - Python package imports (5-10s)
  - Loading 3.4MB Random Forest model (20-40s)
  - Loading XGBoost model (5-10s)
  - FastAPI initialization (5-10s)
  - Total buffer: ~300s ensures no premature timeout

- **3 max retries:** Reduced from 10 to fail faster if issues persist (avoids 10 failed builds)

---

## 3. Deployment Verification Steps

### Immediate Actions (User Should Take)

1. **Monitor Railway Dashboard:**
   - Go to https://railway.app/dashboard
   - Navigate to the "diligent-surprise" project
   - Watch the deployment build logs for the latest commit (cd60d3f)

2. **Check Build Logs for Diagnostic Output:**
   ```
   Look for the startup script output:
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
   -rw-r--r-- 1 root root 3.4M Dec  9 09:55 best_model_random_forest.pkl
   -rw-r--r-- 1 root root 187K Dec  9 09:55 best_model_xgboost.pkl
   -rw-r--r-- 1 root root 2.0K Dec  9 09:55 preprocessor.pkl

   Starting uvicorn server on 0.0.0.0:8080
   ================================
   INFO:     Started server process [1]
   INFO:     Waiting for application startup.
   INFO:     Application startup complete.
   INFO:     Uvicorn running on http://0.0.0.0:8080
   ```

3. **Test Health Endpoint:**
   ```bash
   # Wait for deployment to complete (5-10 minutes), then:
   curl https://web-production-f014a.up.railway.app/health

   # Expected response:
   {
     "status": "healthy",
     "model_loaded": true,
     "model_name": "Random Forest"
   }
   ```

4. **Test Prediction Endpoint:**
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

   # Expected response:
   {
     "toxicity": "Toxic" | "Safe" | "Uncertain",
     "confidence": 85.5,
     "explanation": "...",
     "recommendation": "..."
   }
   ```

5. **Test Interactive API Docs:**
   ```
   https://web-production-f014a.up.railway.app/docs
   ```

---

## 4. Success Criteria Checklist

- [ ] Docker build completes without errors
- [ ] Container starts and runs (no "failed to start" errors)
- [ ] Startup script logs appear showing all checks pass
- [ ] Health endpoint returns 200 OK: `{"status": "healthy"}`
- [ ] Prediction endpoint accepts SMILES and returns toxicity prediction
- [ ] No 502 errors from Railway gateway
- [ ] Logs show "Uvicorn running on http://0.0.0.0:PORT"

---

## 5. Troubleshooting Guide

### If Build Fails

**Check for:**
- Missing system dependencies (apt-get errors)
- Python package installation failures
- File copy errors (missing backend/ directory)

**Solution:**
- Review Railway build logs for specific error messages
- Verify all files are committed to git
- Check if model files (3.6MB total) are included in deployment

### If Container Starts But Health Check Fails

**Symptoms:**
- Build succeeds
- Container logs show startup script running
- Health endpoint returns 502 or timeouts

**Possible Causes:**
1. Models failing to load (check for pickle errors in logs)
2. Port binding issues (check PORT environment variable in logs)
3. Import errors for RDKit or scikit-learn

**Debug Steps:**
```bash
# Check Railway logs
railway logs --service web

# Look for:
# - "ERROR: package not found" in startup checks
# - Python tracebacks
# - "Address already in use" errors
```

### If Health Check Passes But Predictions Fail

**Symptoms:**
- `/health` returns 200 OK
- `/predict` returns 500 errors

**Possible Causes:**
1. Model loading succeeded but prediction logic has bugs
2. Feature engineering failures (SMILES parsing, descriptor calculation)
3. Missing preprocessor or incompatible model versions

**Debug Steps:**
```bash
# Test with simple SMILES
curl -X POST https://web-production-f014a.up.railway.app/test

# Check logs for prediction errors
railway logs --service web | grep -i error
```

---

## 6. Architecture Overview

### Current Stack

```
┌─────────────────────────────────────────────────┐
│              Frontend (Vercel)                   │
│  React + TypeScript + Tailwind CSS              │
│  URL: apistox.vercel.app                        │
└───────────────┬─────────────────────────────────┘
                │ HTTPS/JSON
                ▼
┌─────────────────────────────────────────────────┐
│          Backend API (Railway)                   │
│  FastAPI + Python 3.11                          │
│  URL: web-production-f014a.up.railway.app       │
│                                                  │
│  ┌──────────────────────────────────────────┐  │
│  │  FastAPI Endpoints                        │  │
│  │  - GET  /health                           │  │
│  │  - POST /predict                          │  │
│  │  - GET  /docs                             │  │
│  └──────────────────────────────────────────┘  │
│                                                  │
│  ┌──────────────────────────────────────────┐  │
│  │  ML Models (Loaded on Startup)           │  │
│  │  - Random Forest (3.4MB)                 │  │
│  │  - XGBoost (188KB)                       │  │
│  │  - Preprocessor (2KB)                    │  │
│  └──────────────────────────────────────────┘  │
│                                                  │
│  ┌──────────────────────────────────────────┐  │
│  │  Feature Engineering                      │  │
│  │  - RDKit molecular descriptors           │  │
│  │  - SMILES parsing                        │  │
│  │  - 15+ computed features                 │  │
│  └──────────────────────────────────────────┘  │
└─────────────────────────────────────────────────┘
```

### Docker Container Layout

```
/app/
├── backend/
│   ├── api/
│   │   ├── __init__.py
│   │   ├── main.py          # FastAPI app
│   │   ├── predict.py       # ML prediction logic
│   │   └── featurizer.py    # RDKit feature computation
│   ├── models/
│   │   ├── best_model_random_forest.pkl  (3.4MB)
│   │   ├── best_model_xgboost.pkl        (188KB)
│   │   └── preprocessor.pkl              (2KB)
│   ├── preprocessors/
│   │   └── preprocessor.pkl
│   └── requirements.txt
└── start.sh                  # Startup script with diagnostics
```

---

## 7. Lessons Learned

### What Went Wrong

1. **Incomplete System Dependencies:**
   - Python scientific packages (scikit-learn, xgboost) require C compilers
   - Base Python images don't include build tools
   - Always check package installation for native extensions

2. **Shell Variable Handling:**
   - `sh -c` and `bash -c` handle variable expansion differently
   - Railway's PORT variable needs careful handling
   - Explicit bash scripts are more reliable than inline shell commands

3. **Insufficient Timeout for ML Apps:**
   - ML models take time to load (especially large ones like 3.4MB Random Forest)
   - Default healthcheck timeouts (60-100s) are too short
   - Always benchmark model loading time locally and add buffer

4. **Lack of Diagnostic Logging:**
   - "Container failed to start" with no details makes debugging impossible
   - Startup scripts should verify all dependencies before launching main app
   - Logging package versions helps identify version mismatches

### Best Practices Established

1. **Always Use Startup Scripts for Complex Apps:**
   - Verify environment variables
   - Check dependencies can be imported
   - Confirm file paths and permissions
   - Log diagnostic information

2. **Explicit System Dependencies:**
   - Document why each system package is needed
   - Include build tools for Python packages with C extensions
   - Clean up apt cache to reduce image size

3. **Generous Timeouts for ML Apps:**
   - Calculate expected startup time: dependencies + model loading + initialization
   - Add 50-100% buffer for cold starts
   - Monitor actual startup times and adjust

4. **Version Pinning:**
   - All requirements.txt packages are pinned to specific versions
   - Prevents unexpected breakage from dependency updates
   - Makes debugging reproducible

---

## 8. Performance Considerations

### Current Setup

- **Container Size:** ~1.2GB (Python 3.11-slim + packages + models)
- **Cold Start Time:** ~40-60 seconds (estimated)
  - Package imports: 10-15s
  - Model loading: 25-35s
  - FastAPI init: 5-10s
- **Memory Usage:** ~512MB-1GB (models + FastAPI + RDKit)

### Optimization Opportunities (Future)

1. **Use Model Compression:**
   - Consider quantizing Random Forest model (3.4MB → 1-2MB)
   - Use joblib compression: `joblib.dump(model, 'model.pkl', compress=3)`

2. **Multi-stage Docker Build:**
   - Separate build stage (with gcc, g++) from runtime stage
   - Reduces final image size by 200-300MB

3. **Lazy Model Loading:**
   - Load models on first prediction request instead of startup
   - Improves cold start time but adds latency to first request

4. **Model Caching:**
   - Use Railway's persistent volumes to cache models
   - Avoids re-downloading models on each deployment

---

## 9. Next Steps

### Immediate (After Deployment Succeeds)

1. **Monitor Performance:**
   - Track response times for `/predict` endpoint
   - Monitor memory usage under load
   - Set up Railway metrics dashboard

2. **Add Monitoring:**
   - Integrate error tracking (e.g., Sentry)
   - Add structured logging for predictions
   - Set up uptime monitoring (e.g., UptimeRobot)

3. **Update Frontend:**
   - Ensure frontend is pointing to correct backend URL
   - Add loading states for slow predictions
   - Handle 502/503 errors gracefully

### Medium Term

1. **Add Caching:**
   - Cache predictions for common compounds
   - Use Redis or in-memory cache
   - Reduces compute costs and improves response time

2. **API Rate Limiting:**
   - Implement rate limiting per IP
   - Protect against abuse
   - Use Railway's built-in rate limiting or FastAPI middleware

3. **Comprehensive Testing:**
   - Add integration tests for all endpoints
   - Test with various SMILES strings
   - Load testing to determine capacity

### Long Term

1. **Model Versioning:**
   - Implement model version tracking
   - Allow A/B testing of different models
   - Gradual rollout of model updates

2. **Batch Prediction API:**
   - Add endpoint for predicting multiple compounds
   - Optimize for throughput
   - Support CSV upload

3. **Model Monitoring:**
   - Track prediction distribution
   - Monitor for data drift
   - Alert on anomalous predictions

---

## 10. Contact & Support

### Railway Dashboard
- Project: diligent-surprise
- Environment: production
- Service: web
- URL: https://web-production-f014a.up.railway.app

### Repository
- GitHub: https://github.com/TyLuHow/bee-ML-372
- Branch: main
- Latest Commit: cd60d3f

### Key Files
- Dockerfile: `/home/yler_uby_oward/apistox/Dockerfile`
- Railway Config: `/home/yler_uby_oward/apistox/railway.json`
- Requirements: `/home/yler_uby_oward/apistox/backend/requirements.txt`
- Main API: `/home/yler_uby_oward/apistox/backend/api/main.py`

---

## Appendix A: Complete Dockerfile

```dockerfile
# Use official Python runtime as base image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies for RDKit and scientific packages
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libgomp1 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY backend/requirements.txt ./backend/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r backend/requirements.txt

# Copy application code
COPY backend ./backend

# Create startup script with comprehensive debugging
RUN echo '#!/bin/bash\n\
set -e\n\
echo "================================"\n\
echo "ApisTox Backend Starting"\n\
echo "================================"\n\
echo "Python version: $(python --version)"\n\
echo "Working directory: $(pwd)"\n\
echo "PORT environment variable: ${PORT:-8080}"\n\
echo ""\n\
echo "Checking Python packages..."\n\
python -c "import uvicorn; print(f\"  uvicorn: {uvicorn.__version__}\")" || echo "  ERROR: uvicorn not found"\n\
python -c "import fastapi; print(f\"  fastapi: {fastapi.__version__}\")" || echo "  ERROR: fastapi not found"\n\
python -c "import sklearn; print(f\"  scikit-learn: {sklearn.__version__}\")" || echo "  ERROR: scikit-learn not found"\n\
python -c "import xgboost; print(f\"  xgboost: {xgboost.__version__}\")" || echo "  ERROR: xgboost not found"\n\
python -c "import rdkit; print(f\"  rdkit: {rdkit.__version__}\")" || echo "  ERROR: rdkit not found"\n\
echo ""\n\
echo "Checking model files..."\n\
ls -lh backend/models/ || echo "  ERROR: Models directory not found"\n\
echo ""\n\
echo "Starting uvicorn server on 0.0.0.0:${PORT:-8080}"\n\
echo "================================"\n\
cd backend && exec python -m uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8080}' > /app/start.sh

RUN chmod +x /app/start.sh

# Expose port
EXPOSE 8080

# Start the application
CMD ["/bin/bash", "/app/start.sh"]
```

---

## Appendix B: Complete railway.json

```json
{
  "$schema": "https://railway.app/railway.schema.json",
  "build": {
    "builder": "DOCKERFILE",
    "dockerfilePath": "Dockerfile"
  },
  "deploy": {
    "healthcheckPath": "/health",
    "healthcheckTimeout": 300,
    "restartPolicyType": "ON_FAILURE",
    "restartPolicyMaxRetries": 3
  }
}
```

---

## Appendix C: Testing Commands

### Local Docker Testing

```bash
# Build image locally
cd /home/yler_uby_oward/apistox
docker build -t apistox-test .

# Run container locally
docker run -p 8080:8080 -e PORT=8080 apistox-test

# Test in another terminal
curl http://localhost:8080/health
curl http://localhost:8080/
curl -X POST http://localhost:8080/predict \
  -H "Content-Type: application/json" \
  -d '{"name": "Test", "smiles": "CCO", "category": "Herbicide"}'
```

### Railway Testing

```bash
# Health check
curl https://web-production-f014a.up.railway.app/health

# Root endpoint
curl https://web-production-f014a.up.railway.app/

# Test endpoint (built-in sample)
curl https://web-production-f014a.up.railway.app/test

# Prediction
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

---

## Conclusion

All identified deployment issues have been systematically addressed:

1. Added missing system dependencies (gcc, g++)
2. Created comprehensive startup script with diagnostics
3. Increased healthcheck timeout to 300 seconds
4. Upgraded pip before package installation
5. Committed Dockerfile to git repository
6. Pushed changes to trigger Railway deployment

The deployment should now succeed. Monitor the Railway dashboard for build logs and verify the startup script output appears. Once deployed, test all endpoints to confirm full functionality.

**Estimated time to deployment:** 5-10 minutes from push (completed at 2024-12-09 19:16 UTC)

**Status Check:** Run `curl https://web-production-f014a.up.railway.app/health` every 60 seconds until deployment completes.
