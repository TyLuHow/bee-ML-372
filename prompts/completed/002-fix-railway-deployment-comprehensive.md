<objective>
Comprehensively diagnose and fix ALL Railway deployment issues for the ApisTox backend. The deployment has been failing repeatedly with various errors including missing pip, uvicorn not found, and container startup failures. This prompt will systematically address every issue from the ground up to get a working deployment.
</objective>

<context>
ApisTox is a bee toxicity prediction ML platform with:
- Frontend: React/TypeScript (Vercel)
- Backend: FastAPI + Python ML models (Railway at web-production-f014a.up.railway.app)
- Models: Random Forest (3.4MB) + XGBoost (188KB) + RDKit molecular descriptors

Current situation:
- Multiple failed deployment attempts with different errors
- Build succeeds but container fails to start
- Latest error: "Container failed to start" (no detailed error in logs)
- Reference implementation may exist at /home/yler_uby_oward/tyler-luby-howard-portfolio/Downloads/apistox-pro/

The goal is to achieve a stable, working deployment that serves the ML prediction API.
</context>

<diagnosis_phase>
<title>Comprehensive Root Cause Analysis</title>

Systematically investigate all potential issues:

1. **Examine Reference Implementation**
   - Check if apistox-pro exists and contains working deployment configs
   - Compare Dockerfile, railway.json, requirements.txt with current setup
   - Identify any missing dependencies or configuration differences

2. **Review Current Deployment Configuration**
   - Read @Dockerfile line by line
   - Read @railway.json for deployment settings
   - Read @backend/requirements.txt for dependencies
   - Check if all model files exist in @backend/models/

3. **Analyze Recent Deployment Logs**
   - The build completes successfully (all layers cached)
   - Container fails to start after 28 seconds
   - No clear error message in the logs provided
   - Previous errors included: pip not found, uvicorn not found, "cd" executable not found

4. **Identify Potential Issues**
   - PORT environment variable handling
   - Python module path issues
   - Missing system dependencies for RDKit
   - Incorrect working directory
   - CMD/ENTRYPOINT format problems
   - Health check timing out before app fully loads
</diagnosis_phase>

<fix_strategy>
<title>Systematic Fix Implementation</title>

1. **Create Bulletproof Dockerfile**
   - Use official Python 3.11 base image
   - Install ALL system dependencies needed for RDKit
   - Copy and install requirements properly
   - Use explicit, tested CMD format
   - Add startup logging to debug issues
   - Handle PORT environment variable correctly

2. **Simplify Railway Configuration**
   - Minimal railway.json with just essentials
   - Increase healthcheck timeout for ML model loading
   - Use Dockerfile builder (not Nixpacks)

3. **Verify All Dependencies**
   - Ensure backend/requirements.txt is complete
   - Check model files are committed to git
   - Verify backend/api/__init__.py exists for proper module imports

4. **Add Debugging**
   - Add startup script that logs environment
   - Test that uvicorn can be imported
   - Verify models can be loaded
   - Log PORT variable and binding address

5. **Test Locally First** (if possible)
   - Try building Docker image locally: `docker build -t apistox-test .`
   - Run container locally: `docker run -p 8080:8080 -e PORT=8080 apistox-test`
   - Verify health endpoint responds
</fix_strategy>

<implementation>
Create or update the following files:

1. **Dockerfile** - Rock-solid container configuration:
```dockerfile
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies for RDKit and scientific packages
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libgomp1 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python dependencies
COPY backend/requirements.txt ./backend/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r backend/requirements.txt

# Copy application code
COPY backend ./backend

# Create startup script with debugging
RUN echo '#!/bin/bash\n\
set -e\n\
echo "=== ApisTox Backend Starting ==="\n\
echo "Python version: $(python --version)"\n\
echo "PORT: ${PORT:-8080}"\n\
echo "Working directory: $(pwd)"\n\
echo "Checking uvicorn..."\n\
python -c "import uvicorn; print(f\"uvicorn version: {uvicorn.__version__}\")"\n\
echo "Checking FastAPI..."\n\
python -c "import fastapi; print(f\"fastapi version: {fastapi.__version__}\")"\n\
echo "Checking models..."\n\
ls -lh backend/models/ || echo "Models directory not found"\n\
echo "Starting server on 0.0.0.0:${PORT:-8080}"\n\
cd backend && exec python -m uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8080}' > /app/start.sh

RUN chmod +x /app/start.sh

EXPOSE 8080

CMD ["/bin/bash", "/app/start.sh"]
```

2. **railway.json** - Minimal, robust configuration:
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

3. **Verify backend/api/__init__.py exists** - Create if missing:
```python
"""ApisTox Backend API Package"""
```

4. **Check backend/requirements.txt** - Ensure it includes:
```
fastapi>=0.104.1
uvicorn[standard]>=0.24.0
# ... other dependencies
```
</implementation>

<deployment_steps>
After updating files:

1. Commit changes to git (Railway deploys from git)
2. Deploy to Railway: `railway up --service web`
3. Monitor logs in real-time: `railway logs --service web`
4. Watch for the startup debug messages
5. If it fails, capture the FULL error output
6. Test health endpoint: `curl https://web-production-f014a.up.railway.app/health`
7. If health check fails, increase healthcheckTimeout further
8. Test prediction endpoint with sample SMILES string
</deployment_steps>

<output>
Save detailed findings and fixes to: `./deployment-fix-report.md`

The report should include:
1. **Root Cause Analysis** - What was actually wrong
2. **Changes Made** - File-by-file breakdown of fixes
3. **Deployment Verification** - Proof that it works
4. **Lessons Learned** - How to avoid this in future
5. **Next Steps** - Any remaining improvements needed
</output>

<verification>
Deployment is successful when:
- ✅ Docker build completes without errors
- ✅ Container starts and runs (no "failed to start" errors)
- ✅ Startup script logs appear showing all checks pass
- ✅ Health endpoint returns 200 OK: `{"status": "healthy"}`
- ✅ Prediction endpoint accepts SMILES and returns toxicity prediction
- ✅ No 502 errors from Railway gateway
- ✅ Logs show "Uvicorn running on http://0.0.0.0:PORT"
</verification>

<success_criteria>
- Working Railway deployment serving the ML API
- Health check passing consistently
- Can make successful prediction requests
- Detailed documentation of what was fixed and why
- No more container startup failures
</success_criteria>

<contingency_plans>
If the Dockerfile approach still fails:

**Plan B: Nixpacks with explicit configuration**
- Create working nixpacks.toml based on reference implementation
- Use Railway's Python buildpack with custom start command

**Plan C: Separate build and runtime**
- Use multi-stage Docker build
- Pre-compile Python bytecode
- Minimize container size

**Plan D: Railway template approach**
- Check if Railway has a FastAPI+ML template
- Start from known-working configuration
- Migrate our code incrementally

**Plan E: Alternative platform**
- Document why Railway isn't working
- Recommend Render, Fly.io, or Heroku as alternatives
- Provide migration guide
</contingency_plans>
