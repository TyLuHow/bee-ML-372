# 🎯 ApisTox Command Reference

Quick reference for common commands and operations.

---

## 🚀 First-Time Setup

```bash
# Option 1: Automated (Windows)
QUICK_START.bat

# Option 2: Manual
pip install -r requirements.txt
python src/models.py
python run_comprehensive_analysis.py
```

---

## 🔍 Verification

```bash
# Verify complete setup
python verify_setup.py

# Check Python environment
python --version
pip list | findstr "rdkit xgboost fastapi"

# Check file structure
dir outputs\models
dir outputs\figures
```

---

## 🤖 Model Training

```bash
# Train binary classification model (main)
python src/models.py

# Train ternary classification model
python src/models_ternary.py

# Test model loading
python -c "import joblib; m=joblib.load('outputs/models/best_model.pkl'); print('OK')"
```

---

## 📊 Analysis & Visualization

```bash
# Run all analyses
python run_comprehensive_analysis.py

# Individual analyses
python src/temporal_analysis.py
python src/chemical_space.py
python src/source_comparison.py
python src/toxicophores.py
python src/recommendations.py

# Test scaffold splitting
python test_scaffold_split.py
```

---

## 🌐 API Commands

### Start API

```bash
# Development mode (auto-reload)
python -m uvicorn app.backend.main:app --reload --port 8000

# Production mode
python -m uvicorn app.backend.main:app --host 0.0.0.0 --port 8000

# Different port
python -m uvicorn app.backend.main:app --reload --port 8001
```

### Test API Endpoints

```powershell
# Health check
curl http://localhost:8000/health

# Model info
curl http://localhost:8000/model/info

# Feature importance
curl http://localhost:8000/feature/importance

# Prediction (full features)
curl -X POST http://localhost:8000/predict -H "Content-Type: application/json" -d "@tests/test_prediction.json"

# Prediction (SMILES)
curl -X POST http://localhost:8000/predict/smiles -H "Content-Type: application/json" -d "{\"smiles\":\"CCO\",\"year\":2024,\"herbicide\":0,\"fungicide\":0,\"insecticide\":0,\"other_agrochemical\":1,\"source\":\"PPDB\",\"toxicity_type\":\"Contact\"}"

# Get alternatives
curl http://localhost:8000/recommend/alternatives/CID123

# Toxicophore analysis
curl http://localhost:8000/analysis/toxicophores

# Temporal trends
curl http://localhost:8000/analysis/temporal

# API documentation (in browser)
start http://localhost:8000/docs
```

---

## 🎨 Frontend Commands

```bash
# Navigate to frontend
cd app/frontend

# Install dependencies (first time)
npm install

# Start development server
npm run dev

# Build for production
npm run build

# Preview production build
npm run preview

# Lint code
npm run lint
```

---

## 🧪 Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_api.py -v
pytest tests/test_preprocessing.py -v
pytest tests/test_models.py -v

# Run with coverage
pytest tests/ --cov=src --cov=app.backend --cov-report=term-missing

# Generate HTML coverage report
pytest tests/ --cov=src --cov=app.backend --cov-report=html
start htmlcov/index.html

# Run specific test
pytest tests/test_api.py::TestSMILESEndpoint::test_smiles_prediction_success -v
```

---

## 📦 Data Processing

```bash
# Generate molecular descriptors from SMILES
python src/molecular_features.py

# Run preprocessing pipeline
python src/preprocessing.py

# Validate dataset
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_final.csv'); print(f'Shape: {df.shape}'); print(df.info())"

# Check for missing values
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_final.csv'); print(df.isnull().sum())"
```

---

## 🐛 Debugging

```bash
# Check if port 8000 is in use
netstat -ano | findstr :8000

# Kill process on port 8000 (replace PID)
taskkill /PID <PID> /F

# Check logs (if API crashes)
python -m uvicorn app.backend.main:app --log-level debug

# Test model prediction directly
python -c "import joblib; import numpy as np; m=joblib.load('outputs/models/best_model.pkl'); print(m.predict([[1]*15]))"

# Verify RDKit SMILES parsing
python -c "from rdkit import Chem; mol=Chem.MolFromSmiles('CCO'); print('Valid' if mol else 'Invalid')"
```

---

## 📝 Git Commands

```bash
# Check status
git status

# Add changes
git add .

# Commit
git commit -m "Description of changes"

# Push to remote
git push origin master

# Pull latest changes
git pull origin master

# View commit history
git log --oneline -10
```

---

## 🗂️ File Management

```bash
# List all output files
dir /S outputs

# Check output file sizes
dir outputs\models
dir outputs\figures
dir outputs\analysis

# Clean outputs (be careful!)
# rmdir /S /Q outputs\figures
# rmdir /S /Q outputs\analysis

# Create missing directories
mkdir outputs\models
mkdir outputs\figures
mkdir outputs\analysis
```

---

## 🚢 Deployment

```bash
# Build Docker image (if using Docker)
docker build -t apistox .

# Run Docker container
docker run -p 8000:8000 apistox

# Deploy to Vercel (from project root)
vercel --prod

# Check deployment logs
vercel logs
```

---

## 📊 Generate Reports

```bash
# Generate all reports
python src/toxicophores.py          # → outputs/TOXICOPHORE_REPORT.md
python src/recommendations.py       # → outputs/ALTERNATIVES_REPORT.md
python run_comprehensive_analysis.py # → outputs/TEMPORAL_ANALYSIS_REPORT.md

# View reports
type outputs\TOXICOPHORE_REPORT.md
type outputs\ALTERNATIVES_REPORT.md
type outputs\TEMPORAL_ANALYSIS_REPORT.md
```

---

## 🔧 Maintenance

```bash
# Update dependencies
pip install --upgrade -r requirements.txt

# Check for outdated packages
pip list --outdated

# Reinstall specific package
pip install --force-reinstall rdkit==2023.9.5

# Clear Python cache
del /S /Q __pycache__
del /S /Q *.pyc

# Retrain all models
python src/models.py
python src/models_ternary.py
```

---

## 💡 Tips & Tricks

### Run Multiple Commands

```batch
REM Create a batch file for common workflow
@echo off
python src/models.py && python run_comprehensive_analysis.py && python -m uvicorn app.backend.main:app --reload
```

### Set Environment Variables

```powershell
# PowerShell
$env:API_URL = "http://localhost:8000"
$env:DEBUG = "True"

# Command Prompt
set API_URL=http://localhost:8000
set DEBUG=True
```

### Monitor API Logs

```bash
# Run API with detailed logging
python -m uvicorn app.backend.main:app --reload --log-level debug > api.log 2>&1

# Watch log file (PowerShell)
Get-Content api.log -Wait
```

---

## 🆘 Emergency Commands

```bash
# Reset everything (DANGER!)
# rmdir /S /Q outputs
# rmdir /S /Q __pycache__
# pip uninstall -y -r requirements.txt

# Start fresh
git clean -fdx  # Remove all untracked files
pip install -r requirements.txt
python QUICK_START.bat
```

---

## 📱 Shortcuts

```bash
# Alias for common commands (PowerShell profile)
function Start-API { python -m uvicorn app.backend.main:app --reload --port 8000 }
function Start-Frontend { cd app/frontend; npm run dev }
function Run-Tests { pytest tests/ -v }

# Usage
Start-API
Start-Frontend
Run-Tests
```

---

## 🎓 Learning Resources

```bash
# View API documentation
start http://localhost:8000/docs

# Check model performance
python -c "import joblib; m=joblib.load('outputs/models/best_model.pkl'); print(m)"

# Explore dataset
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_final.csv'); print(df.describe())"

# View analysis results
start outputs\figures\temporal_trends.png
start outputs\TOXICOPHORE_REPORT.md
```

---

*For detailed instructions, see `SETUP_GUIDE.md`*  
*For API documentation, visit http://localhost:8000/docs*




