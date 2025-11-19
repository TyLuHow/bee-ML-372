# ApisTox Setup & Deployment Guide

**Document Version**: 1.0
**Last Updated**: November 19, 2025
**Platforms**: Linux, macOS, Windows (WSL)

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Local Development Setup](#local-development-setup)
3. [Running the Application](#running-the-application)
4. [Docker Deployment](#docker-deployment)
5. [Cloud Deployment](#cloud-deployment)
6. [Serverless Deployment](#serverless-deployment)
7. [Production Configuration](#production-configuration)
8. [Troubleshooting](#troubleshooting)
9. [Maintenance](#maintenance)

---

## Prerequisites

### System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **CPU** | 2 cores | 4+ cores |
| **RAM** | 4 GB | 8+ GB |
| **Disk** | 2 GB free | 5+ GB free |
| **OS** | Ubuntu 20.04, macOS 11, Windows 10 (WSL2) | Latest LTS |

### Software Dependencies

#### Backend Requirements

```bash
# Python 3.9 or higher
python3 --version  # Should show 3.9.x or 3.10.x or 3.11.x

# pip package manager
pip --version

# Virtual environment (optional but recommended)
python3 -m venv --help
```

#### Frontend Requirements

```bash
# Node.js 18 or higher
node --version  # Should show v18.x.x or v20.x.x

# npm package manager
npm --version  # Should show 9.x.x or higher
```

#### Optional (for Docker deployment)

```bash
# Docker
docker --version  # Should show 20.x.x or higher

# Docker Compose
docker-compose --version  # Should show 2.x.x or higher
```

---

## Local Development Setup

### Step 1: Clone Repository

```bash
# Clone the repository
git clone https://github.com/YOUR_ORG/apistox.git
cd apistox

# Verify directory structure
ls -la
# Expected: app/, src/, data/, outputs/, tests/, docs/, etc.
```

### Step 2: Backend Setup

#### 2.1 Create Python Virtual Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
# On Linux/macOS:
source venv/bin/activate
# On Windows (WSL):
source venv/bin/activate
# On Windows (CMD):
venv\Scripts\activate.bat
# On Windows (PowerShell):
venv\Scripts\Activate.ps1

# Verify activation (should show (venv) in prompt)
which python  # Should point to venv/bin/python
```

#### 2.2 Install Python Dependencies

```bash
# Upgrade pip
pip install --upgrade pip

# Install all dependencies (development + production)
pip install -r requirements.txt

# This installs:
# - FastAPI, Uvicorn (API server)
# - scikit-learn, XGBoost, LightGBM (ML)
# - RDKit (cheminformatics)
# - Pandas, NumPy (data processing)
# - SHAP, LIME (interpretability)
# - Matplotlib, Seaborn, Plotly (visualization)
# - imbalanced-learn (SMOTE)
# - pytest (testing)

# Verify installation
python -c "import fastapi, xgboost, rdkit; print('✓ All imports successful')"
```

**Troubleshooting RDKit Installation**:

If RDKit fails to install via pip, use conda:

```bash
# Install conda (if not already installed)
# Download from: https://docs.conda.io/en/latest/miniconda.html

# Create conda environment
conda create -n apistox python=3.9
conda activate apistox

# Install RDKit via conda
conda install -c conda-forge rdkit

# Install other dependencies via pip
pip install -r requirements.txt
```

#### 2.3 Verify Data and Models

```bash
# Check if dataset exists
ls -lh data/raw/dataset_with_descriptors.csv
# Expected: ~0.8 MB file, 1035 rows

# Check if models exist (if provided)
ls -lh outputs/models/
# Expected: best_model_xgboost.pkl (~13 MB), preprocessor in outputs/preprocessors/

# If models don't exist, train them:
python src/preprocessing.py  # Creates preprocessor.pkl (~30 seconds)
python src/models.py         # Creates all models (~5 minutes)
```

#### 2.4 Run Backend Tests (Optional)

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_api.py -v

# Run with coverage report
pytest tests/ --cov=src --cov=app.backend --cov-report=html
# Open htmlcov/index.html in browser to see coverage
```

### Step 3: Frontend Setup

#### 3.1 Install Node.js Dependencies

```bash
# Navigate to frontend directory
cd app/frontend

# Install dependencies (takes 2-3 minutes)
npm install

# This installs:
# - react, react-dom (UI framework)
# - typescript (type safety)
# - vite (build tool)
# - tailwindcss (styling)
# - axios (HTTP client)
# - recharts (charts)

# Verify installation
npm list --depth=0
```

#### 3.2 Configure API URL

```bash
# Create environment file (optional, for production)
echo "VITE_API_URL=http://localhost:8000" > .env

# For development, default is http://localhost:8000 (no .env needed)
```

#### 3.3 Build Frontend (Optional, for testing)

```bash
# Build production bundle
npm run build

# Output will be in dist/ directory
ls -lh dist/
# Expected: index.html, assets/ (JS/CSS bundles)

# Preview production build
npm run preview
# Opens at http://localhost:4173
```

---

## Running the Application

### Option 1: Development Mode (Hot Reload)

#### Terminal 1: Backend Server

```bash
# From project root, with venv activated
cd /path/to/apistox
source venv/bin/activate  # Or conda activate apistox

# Run FastAPI with auto-reload
uvicorn app.backend.main:app --reload --port 8000

# Alternative: Run directly with Python
python app/backend/main.py

# Expected output:
# INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
# INFO:     Started reloader process [12345] using WatchFiles
# INFO:     Started server process [12346]
# INFO:     Waiting for application startup.
# INFO:     Application startup complete.

# Test API is running:
curl http://localhost:8000/health
# Expected: {"status":"healthy","model_loaded":true,"timestamp":"..."}
```

**API Endpoints Available**:
- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc
- **OpenAPI JSON**: http://localhost:8000/openapi.json

#### Terminal 2: Frontend Dev Server

```bash
# From project root
cd app/frontend

# Run Vite dev server (port 5173)
npm run dev

# Expected output:
# VITE v5.0.8  ready in 523 ms
#
# ➜  Local:   http://localhost:5173/
# ➜  Network: use --host to expose
# ➜  press h to show help

# Open browser to: http://localhost:5173
```

**Features in Dev Mode**:
- ✅ Hot Module Reload (HMR) for frontend
- ✅ Auto-reload for backend code changes
- ✅ Source maps for debugging
- ✅ Fast refresh (<100ms)

### Option 2: Production Mode (Single Process)

```bash
# Build frontend
cd app/frontend
npm run build
cd ../..

# Serve frontend from backend (if configured)
# OR use separate web server (Nginx, Apache)

# Run backend with production settings
gunicorn app.backend.main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000

# Explanation:
# -w 4: 4 worker processes (adjust to CPU cores)
# -k uvicorn.workers.UvicornWorker: ASGI worker class
# --bind 0.0.0.0:8000: Listen on all interfaces
```

### Option 3: Using Screen/Tmux (Background Processes)

```bash
# Start backend in background
screen -S apistox-backend
source venv/bin/activate
uvicorn app.backend.main:app --host 0.0.0.0 --port 8000
# Press Ctrl+A, then D to detach

# Start frontend in background
screen -S apistox-frontend
cd app/frontend
npm run dev -- --host 0.0.0.0
# Press Ctrl+A, then D to detach

# List running screens
screen -ls

# Reattach to screen
screen -r apistox-backend
```

---

## Docker Deployment

### Option 1: Docker Compose (Recommended)

#### Step 1: Build and Run

```bash
# From project root
docker-compose up --build

# OR run in detached mode (background)
docker-compose up -d --build

# Expected output:
# Building backend...
# Building frontend...
# Creating apistox_backend_1 ... done
# Creating apistox_frontend_1 ... done

# View logs
docker-compose logs -f

# View specific service logs
docker-compose logs -f backend
docker-compose logs -f frontend
```

#### Step 2: Access Application

```
Frontend: http://localhost:80
Backend API: http://localhost:8000
Swagger Docs: http://localhost:8000/docs
```

#### Step 3: Stop and Clean Up

```bash
# Stop containers
docker-compose down

# Stop and remove volumes
docker-compose down -v

# Stop and remove images
docker-compose down --rmi all
```

### Option 2: Individual Docker Containers

#### Backend Container

```bash
# Build backend image
docker build -f Dockerfile.backend -t apistox-backend:latest .

# Run backend container
docker run -d \
  --name apistox-backend \
  -p 8000:8000 \
  -v $(pwd)/outputs:/app/outputs \
  apistox-backend:latest

# View logs
docker logs -f apistox-backend

# Stop and remove
docker stop apistox-backend
docker rm apistox-backend
```

#### Frontend Container

```bash
# Build frontend image
docker build -f Dockerfile.frontend -t apistox-frontend:latest .

# Run frontend container
docker run -d \
  --name apistox-frontend \
  -p 80:80 \
  apistox-frontend:latest

# View logs
docker logs -f apistox-frontend

# Stop and remove
docker stop apistox-frontend
docker rm apistox-frontend
```

### Docker Best Practices

1. **Multi-stage builds**: Already implemented in Dockerfiles
2. **Layer caching**: Order commands from least to most frequently changed
3. **Secrets**: Use Docker secrets or environment variables (never commit secrets)
4. **Health checks**: Add to docker-compose.yml

```yaml
# Add to docker-compose.yml
services:
  backend:
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 40s
```

---

## Cloud Deployment

### AWS Deployment

#### Option 1: EC2 Instance

**Step 1: Launch EC2 Instance**

```bash
# Instance type: t3.medium or larger
# OS: Ubuntu 22.04 LTS
# Storage: 20 GB SSD
# Security Group: Allow ports 22 (SSH), 80 (HTTP), 443 (HTTPS), 8000 (API)
```

**Step 2: Connect and Setup**

```bash
# SSH into instance
ssh -i your-key.pem ubuntu@ec2-xx-xx-xx-xx.compute.amazonaws.com

# Update system
sudo apt update && sudo apt upgrade -y

# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker ubuntu
# Log out and log back in

# Install Docker Compose
sudo curl -L "https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# Clone repository
git clone https://github.com/YOUR_ORG/apistox.git
cd apistox

# Run with Docker Compose
docker-compose up -d
```

**Step 3: Configure Nginx (Optional, for reverse proxy)**

```bash
# Install Nginx
sudo apt install nginx -y

# Create Nginx config
sudo nano /etc/nginx/sites-available/apistox

# Add configuration:
server {
    listen 80;
    server_name your-domain.com;

    # Frontend
    location / {
        proxy_pass http://localhost:80;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    # Backend API
    location /api {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}

# Enable site
sudo ln -s /etc/nginx/sites-available/apistox /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

**Step 4: SSL Certificate (Let's Encrypt)**

```bash
# Install Certbot
sudo apt install certbot python3-certbot-nginx -y

# Get certificate
sudo certbot --nginx -d your-domain.com

# Auto-renewal is configured automatically
```

#### Option 2: ECS (Elastic Container Service)

```bash
# 1. Push images to ECR
aws ecr create-repository --repository-name apistox-backend
aws ecr create-repository --repository-name apistox-frontend

# 2. Build and push
docker build -f Dockerfile.backend -t apistox-backend .
docker tag apistox-backend:latest <account-id>.dkr.ecr.<region>.amazonaws.com/apistox-backend:latest
docker push <account-id>.dkr.ecr.<region>.amazonaws.com/apistox-backend:latest

# 3. Create ECS task definition (use AWS Console or CLI)
# 4. Create ECS service with load balancer
# 5. Configure auto-scaling
```

#### Option 3: Lambda + API Gateway (Serverless)

- See "Serverless Deployment" section below

### Google Cloud Platform (GCP)

#### Option 1: Cloud Run (Recommended)

```bash
# 1. Install gcloud CLI
curl https://sdk.cloud.google.com | bash
gcloud init

# 2. Build and push to Container Registry
gcloud builds submit --tag gcr.io/PROJECT_ID/apistox-backend
gcloud builds submit --tag gcr.io/PROJECT_ID/apistox-frontend

# 3. Deploy to Cloud Run
gcloud run deploy apistox-backend \
  --image gcr.io/PROJECT_ID/apistox-backend \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2

gcloud run deploy apistox-frontend \
  --image gcr.io/PROJECT_ID/apistox-frontend \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated

# 4. Get URLs
gcloud run services list
```

#### Option 2: Compute Engine (VM)

- Similar to AWS EC2 steps above

### Microsoft Azure

#### Option 1: App Service

```bash
# 1. Install Azure CLI
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash

# 2. Login
az login

# 3. Create resource group
az group create --name apistox-rg --location eastus

# 4. Create App Service plan
az appservice plan create --name apistox-plan --resource-group apistox-rg --sku B1 --is-linux

# 5. Deploy backend
az webapp create --resource-group apistox-rg --plan apistox-plan --name apistox-backend --deployment-container-image-name apistox-backend:latest

# 6. Deploy frontend
az webapp create --resource-group apistox-rg --plan apistox-plan --name apistox-frontend --deployment-container-image-name apistox-frontend:latest
```

---

## Serverless Deployment

### Vercel (Frontend + Serverless Functions)

#### Step 1: Install Vercel CLI

```bash
npm install -g vercel
```

#### Step 2: Configure Project

```bash
# From project root
vercel login

# Initialize project
vercel

# Follow prompts:
# - Link to existing project? N (first time)
# - Project name: apistox
# - Directory: ./
# - Override settings? N
```

#### Step 3: Configure vercel.json (Already Present)

```json
{
  "version": 2,
  "builds": [
    {
      "src": "app/backend/main.py",
      "use": "@vercel/python",
      "config": {
        "maxLambdaSize": "15mb"
      }
    },
    {
      "src": "app/frontend/package.json",
      "use": "@vercel/static-build",
      "config": {
        "distDir": "dist"
      }
    }
  ],
  "routes": [
    {
      "src": "/api/(.*)",
      "dest": "app/backend/main.py"
    },
    {
      "src": "/(.*)",
      "dest": "app/frontend/$1"
    }
  ]
}
```

#### Step 4: Deploy

```bash
# Deploy to production
vercel --prod

# Expected output:
# 🔍  Inspect: https://vercel.com/...
# ✅  Production: https://apistox.vercel.app

# View logs
vercel logs
```

**Limitations**:
- 15 MB function size limit (exclude SHAP)
- 10 second execution timeout (sufficient for predictions)
- Stateless (use cloud storage for models)

### Railway

#### Step 1: Install Railway CLI

```bash
npm install -g @railway/cli
```

#### Step 2: Deploy

```bash
# Login
railway login

# Initialize project
railway init

# Deploy
railway up

# Get URL
railway open
```

### Heroku

#### Step 1: Install Heroku CLI

```bash
curl https://cli-assets.heroku.com/install.sh | sh
```

#### Step 2: Deploy

```bash
# Login
heroku login

# Create app
heroku create apistox

# Add Procfile (already present)
# web: gunicorn app.backend.main:app -w 4 -k uvicorn.workers.UvicornWorker

# Deploy
git push heroku main

# Open app
heroku open
```

---

## Production Configuration

### Environment Variables

Create `.env` file (never commit to git):

```bash
# Backend
API_HOST=0.0.0.0
API_PORT=8000
DEBUG=false
ALLOWED_ORIGINS=https://your-frontend-domain.com

# Database (future)
DATABASE_URL=postgresql://user:password@localhost:5432/apistox

# Redis (future)
REDIS_URL=redis://localhost:6379

# Secrets
SECRET_KEY=your-secret-key-here
API_KEY_SALT=your-salt-here

# Model paths
MODEL_PATH=/app/outputs/models/best_model_xgboost.pkl
PREPROCESSOR_PATH=/app/outputs/preprocessors/preprocessor.pkl

# Monitoring
SENTRY_DSN=https://...@sentry.io/...

# Frontend
VITE_API_URL=https://api.your-domain.com
```

### Load Environment Variables

**Backend (Python)**:

```python
# Install python-dotenv
pip install python-dotenv

# In main.py
from dotenv import load_dotenv
import os

load_dotenv()

API_HOST = os.getenv("API_HOST", "127.0.0.1")
API_PORT = int(os.getenv("API_PORT", 8000))
DEBUG = os.getenv("DEBUG", "false").lower() == "true"
```

**Frontend (Vite)**:

Vite automatically loads `.env` files. Access with `import.meta.env.VITE_*`:

```typescript
// In api.ts
const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000';
```

### Security Hardening

1. **HTTPS Only** (use Nginx with SSL):

```nginx
server {
    listen 443 ssl http2;
    ssl_certificate /etc/letsencrypt/live/your-domain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/your-domain.com/privkey.pem;

    # Redirect HTTP to HTTPS
    if ($scheme != "https") {
        return 301 https://$server_name$request_uri;
    }
}
```

2. **CORS Configuration** (backend):

```python
# In main.py
app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://your-frontend-domain.com"],  # Specific origins only
    allow_credentials=True,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)
```

3. **Rate Limiting** (add to backend):

```bash
pip install slowapi

# In main.py
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

limiter = Limiter(key_func=get_remote_address)
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

@app.post("/predict")
@limiter.limit("100/minute")
async def predict(request: Request, input_data: PredictionInput):
    ...
```

4. **Authentication** (add API keys):

```python
from fastapi.security import APIKeyHeader

API_KEY_HEADER = APIKeyHeader(name="X-API-Key")

def verify_api_key(api_key: str = Depends(API_KEY_HEADER)):
    if api_key != os.getenv("API_KEY"):
        raise HTTPException(status_code=403, detail="Invalid API key")
    return api_key

@app.post("/predict")
async def predict(input_data: PredictionInput, api_key: str = Depends(verify_api_key)):
    ...
```

### Database Setup (Future Enhancement)

```bash
# Install PostgreSQL
sudo apt install postgresql postgresql-contrib

# Create database
sudo -u postgres createdb apistox

# Create user
sudo -u postgres createuser apistox_user

# Grant privileges
sudo -u postgres psql
ALTER USER apistox_user WITH PASSWORD 'your-password';
GRANT ALL PRIVILEGES ON DATABASE apistox TO apistox_user;

# Install Python client
pip install psycopg2-binary sqlalchemy

# Create tables (use Alembic for migrations)
pip install alembic
alembic init migrations
```

---

## Troubleshooting

### Common Issues

#### 1. Port Already in Use

**Problem**: `uvicorn.error ERROR: [Errno 48] Address already in use`

**Solution**:
```bash
# Find process using port 8000
lsof -i :8000
# OR
netstat -vanp tcp | grep 8000

# Kill process
kill -9 <PID>

# Use different port
uvicorn app.backend.main:app --port 8001
```

#### 2. ModuleNotFoundError: No module named 'rdkit'

**Problem**: RDKit not installed

**Solution**:
```bash
# If using pip failed, use conda
conda install -c conda-forge rdkit

# OR use pre-built wheel (Linux)
pip install rdkit-pypi
```

#### 3. Model File Not Found

**Problem**: `FileNotFoundError: [Errno 2] No such file or directory: 'outputs/models/best_model_xgboost.pkl'`

**Solution**:
```bash
# Train models
python src/preprocessing.py
python src/models.py

# Verify files created
ls -lh outputs/models/
ls -lh outputs/preprocessors/
```

#### 4. CORS Error in Browser

**Problem**: `Access to XMLHttpRequest... has been blocked by CORS policy`

**Solution**:
```python
# In app/backend/main.py, update CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:3000"],  # Add your frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

#### 5. npm install Fails

**Problem**: `EACCES: permission denied`

**Solution**:
```bash
# Fix npm permissions
mkdir ~/.npm-global
npm config set prefix '~/.npm-global'
echo 'export PATH=~/.npm-global/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

# Retry
npm install
```

#### 6. Docker Build Fails (Out of Memory)

**Problem**: Docker build killed due to OOM

**Solution**:
```bash
# Increase Docker memory limit (Docker Desktop)
# Settings → Resources → Memory → 4 GB or higher

# OR build with reduced parallelism
docker-compose build --parallel 1
```

### Logging and Debugging

#### Backend Logs

```bash
# Development (console output)
uvicorn app.backend.main:app --reload --log-level debug

# Production (log to file)
uvicorn app.backend.main:app --log-config logging.conf

# Docker logs
docker-compose logs -f backend
```

#### Frontend Logs

```bash
# Browser console (F12)
# Check for:
# - Network errors (failed API calls)
# - JavaScript errors
# - CORS issues

# Vite dev server logs
npm run dev
```

#### Model Predictions Logging

```python
# Add to app/backend/main.py
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@app.post("/predict")
async def predict(input_data: PredictionInput):
    logger.info(f"Received prediction request: {input_data.dict()}")
    # ... prediction logic
    logger.info(f"Prediction result: {result}")
    return result
```

---

## Maintenance

### Regular Tasks

#### Daily
- ✅ Check application health (`/health` endpoint)
- ✅ Monitor error logs
- ✅ Check disk space (`df -h`)

#### Weekly
- ✅ Review prediction history
- ✅ Update dependencies (`pip list --outdated`, `npm outdated`)
- ✅ Backup database (future)
- ✅ Check SSL certificate expiry

#### Monthly
- ✅ Update Docker images
- ✅ Security patches (`apt update && apt upgrade`)
- ✅ Review and rotate API keys
- ✅ Test backup restoration

#### Quarterly
- ✅ Retrain models with new data
- ✅ Performance benchmarking
- ✅ Capacity planning
- ✅ Dependency major version upgrades

### Backup Strategy

```bash
# Backup script (save as backup.sh)
#!/bin/bash
BACKUP_DIR="/backups/apistox/$(date +%Y%m%d)"
mkdir -p $BACKUP_DIR

# Backup data
cp data/raw/*.csv $BACKUP_DIR/

# Backup models
cp -r outputs/models/ $BACKUP_DIR/
cp -r outputs/preprocessors/ $BACKUP_DIR/

# Backup database (future)
# pg_dump apistox > $BACKUP_DIR/database.sql

# Compress
tar -czf $BACKUP_DIR.tar.gz $BACKUP_DIR
rm -rf $BACKUP_DIR

# Upload to cloud storage
# aws s3 cp $BACKUP_DIR.tar.gz s3://apistox-backups/
```

### Monitoring

```bash
# Install monitoring tools
pip install prometheus-fastapi-instrumentator

# Add to main.py
from prometheus_fastapi_instrumentator import Instrumentator

app = FastAPI()
Instrumentator().instrument(app).expose(app)

# Metrics available at: http://localhost:8000/metrics
```

### Scaling

**Horizontal Scaling** (add more instances):
```bash
# Docker Compose
docker-compose up --scale backend=3

# Kubernetes (future)
kubectl scale deployment apistox-backend --replicas=3
```

**Vertical Scaling** (larger instances):
- Increase CPU/RAM in cloud provider console
- Update docker-compose.yml resource limits

---

## Quick Reference

### Essential Commands

| Task | Command |
|------|---------|
| Start backend (dev) | `uvicorn app.backend.main:app --reload` |
| Start frontend (dev) | `cd app/frontend && npm run dev` |
| Run tests | `pytest tests/ -v` |
| Build frontend | `cd app/frontend && npm run build` |
| Docker up | `docker-compose up -d` |
| Docker down | `docker-compose down` |
| View logs | `docker-compose logs -f` |
| Train models | `python src/models.py` |
| Check health | `curl http://localhost:8000/health` |

### Useful URLs

| Service | URL |
|---------|-----|
| Frontend (dev) | http://localhost:5173 |
| Backend API | http://localhost:8000 |
| API Docs (Swagger) | http://localhost:8000/docs |
| API Docs (ReDoc) | http://localhost:8000/redoc |
| Metrics (future) | http://localhost:8000/metrics |

---

**Document Maintained By**: ApisTox DevOps Team
**Last Review**: November 19, 2025
**Next Review**: December 2025
