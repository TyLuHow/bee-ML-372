# 🐝 ApisTox Complete Setup Guide

**Last Updated**: November 8, 2025  
**Estimated Time**: 30-45 minutes  
**Skill Level**: Intermediate

---

## 📋 Prerequisites

- **Python 3.9-3.12** installed
- **Node.js 16+** and npm installed
- **Git** installed
- **8GB RAM minimum** (16GB recommended)
- **5GB free disk space**

---

## 🚀 Quick Start (5 Steps)

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Run data processing & model training
python src/preprocessing.py
python src/models.py

# 3. Run all analyses
python run_comprehensive_analysis.py

# 4. Start the API (in one terminal)
python -m uvicorn app.backend.main:app --reload --port 8000

# 5. Start the frontend (in another terminal)
cd app/frontend
npm install
npm run dev
```

---

## 📖 Detailed Setup Instructions

### Step 1: Environment Setup

#### Option A: Using pip (Recommended)

```bash
# Navigate to project directory
cd C:\Users\Tyler Luby Howard\apis_tox_dataset

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "import rdkit; import xgboost; import fastapi; print('✓ All packages installed')"
```

#### Option B: Using conda (Alternative)

```bash
# Create conda environment
conda create -n apistox python=3.10 -y
conda activate apistox

# Install packages
pip install -r requirements.txt
```

---

### Step 2: Data Processing & Feature Engineering

#### 2.1 Generate Molecular Descriptors

```bash
# This adds RDKit molecular descriptors to the dataset
python src/molecular_features.py

# Expected output:
# ✓ Processed 1035 compounds
# ✓ Generated 15 molecular descriptors
# ✓ Saved to data/raw/dataset_with_descriptors.csv
```

#### 2.2 Run Preprocessing

```bash
# Prepare features for modeling
python src/preprocessing.py

# Expected output:
# ✓ Features prepared
# ✓ Train/val/test split created
# ✓ Preprocessor saved to outputs/models/preprocessor.pkl
```

**Verify**: Check that `outputs/models/preprocessor.pkl` exists

---

### Step 3: Model Training

#### 3.1 Train Binary Classification Model (EPA Label)

```bash
python src/models.py

# This trains multiple models and selects the best one
# Expected time: 5-10 minutes
# Expected output:
# ✓ Trained 5 models
# ✓ Best model: XGBoost (F1: 0.84)
# ✓ Saved to outputs/models/best_model.pkl
```

#### 3.2 Train Ternary Classification Model (PPDB Levels) - OPTIONAL

```bash
python src/models_ternary.py

# Expected output:
# ✓ Best model saved to outputs/models/best_model_ternary_xgboost.pkl
# ✓ Confusion matrix saved
```

**Verify**: Check that `outputs/models/best_model.pkl` exists

---

### Step 4: Run All Analyses

#### 4.1 Comprehensive Analysis (Temporal + Chemical Space)

```bash
python run_comprehensive_analysis.py

# Expected time: 2-3 minutes
# Generates:
# - Temporal trend analysis
# - Chemical space visualizations (PCA, t-SNE)
# - 10+ figures in outputs/figures/
```

#### 4.2 Source Comparison (ECOTOX vs PPDB)

```bash
python src/source_comparison.py

# Expected output:
# ✓ Chi-square test completed
# ✓ Source comparison plots saved
# ✓ outputs/figures/source_comparison.png
```

#### 4.3 Toxicophore Identification

```bash
python src/toxicophores.py

# Expected time: 3-5 minutes
# Expected output:
# ✓ Identified 20 toxicophores
# ✓ Statistical analysis complete
# ✓ outputs/TOXICOPHORE_REPORT.md created
```

#### 4.4 Alternative Compound Recommendations

```bash
python src/recommendations.py

# Expected time: 2-3 minutes
# Expected output:
# ✓ Generated alternatives for all toxic compounds
# ✓ outputs/ALTERNATIVES_REPORT.md created
# ✓ outputs/analysis/alternatives.csv
```

**Verify**: Check `outputs/figures/` has 20+ PNG files

---

### Step 5: Start the Backend API

#### 5.1 Test Model Loading

```bash
# Quick test to ensure models load correctly
python -c "import joblib; model = joblib.load('outputs/models/best_model.pkl'); print('✓ Model loaded successfully')"
```

#### 5.2 Start API Server

```bash
# Start FastAPI server
python -m uvicorn app.backend.main:app --reload --port 8000

# Expected output:
# ✓ Model loaded from outputs/models/best_model.pkl
# ✓ Preprocessor loaded from outputs/models/preprocessor.pkl
# INFO: Uvicorn running on http://127.0.0.1:8000
```

**Keep this terminal open!**

#### 5.3 Test API Endpoints (in a new terminal)

```powershell
# Test health endpoint
curl http://localhost:8000/health

# Test prediction endpoint (use the test file)
curl -X POST http://localhost:8000/predict -H "Content-Type: application/json" -d "@tests/test_prediction.json"

# Test SMILES endpoint
curl -X POST http://localhost:8000/predict/smiles -H "Content-Type: application/json" -d '{\"smiles\":\"CCO\",\"year\":2024,\"herbicide\":0,\"fungicide\":0,\"insecticide\":0,\"other_agrochemical\":1,\"source\":\"PPDB\",\"toxicity_type\":\"Contact\"}'
```

**Expected**: All endpoints return 200 OK with JSON responses

---

### Step 6: Start the Frontend

#### 6.1 Install Frontend Dependencies

```bash
# Navigate to frontend directory
cd app/frontend

# Install packages (first time only)
npm install

# Expected time: 2-3 minutes
# Expected output:
# added 200+ packages
```

#### 6.2 Configure API URL (if needed)

```bash
# Check if .env file exists
type .env

# If not, create it:
echo REACT_APP_API_URL=http://localhost:8000 > .env
```

#### 6.3 Start Development Server

```bash
# Still in app/frontend directory
npm run dev

# Expected output:
# VITE v5.x.x ready in 500 ms
# ➜ Local: http://localhost:5173/
# ➜ Network: use --host to expose
```

**Keep this terminal open!**

#### 6.4 Test Frontend

Open your browser to **http://localhost:5173**

You should see:
- 🐝 ApisTox header
- Prediction form with inputs
- Model info panel

Try making a prediction!

---

### Step 7: Run Tests

#### 7.1 Unit Tests

```bash
# Run all tests
pytest tests/ -v

# Expected output:
# test_preprocessing.py::test_data_cleaning PASSED
# test_models.py::test_model_training PASSED
# test_api.py::test_health_check PASSED
# ... 15+ tests PASSED
```

#### 7.2 API Integration Tests

```bash
# Run API-specific tests
pytest tests/test_api.py -v

# Expected output:
# test_health_check PASSED
# test_smiles_prediction_success PASSED
# test_prediction_success PASSED
# ... 8+ tests PASSED
```

#### 7.3 Coverage Report

```bash
pytest tests/ --cov=src --cov=app.backend --cov-report=html

# Opens coverage report in browser
start htmlcov/index.html
```

---

## 🎯 Verification Checklist

After completing all steps, verify:

### Files Created

- [ ] `outputs/models/best_model.pkl` (trained model)
- [ ] `outputs/models/preprocessor.pkl` (preprocessor)
- [ ] `outputs/figures/` contains 20+ PNG plots
- [ ] `outputs/analysis/` contains JSON result files
- [ ] `outputs/TOXICOPHORE_REPORT.md` exists
- [ ] `outputs/ALTERNATIVES_REPORT.md` exists
- [ ] `outputs/TEMPORAL_ANALYSIS_REPORT.md` exists

### Services Running

- [ ] Backend API at http://localhost:8000
- [ ] API docs at http://localhost:8000/docs
- [ ] Frontend at http://localhost:5173
- [ ] Can make predictions through UI

### Tests Passing

- [ ] `pytest tests/` shows 90%+ pass rate
- [ ] API endpoints respond correctly
- [ ] Frontend connects to backend

---

## 🐛 Troubleshooting

### Issue 1: "ModuleNotFoundError: No module named 'rdkit'"

**Solution**:
```bash
pip install rdkit==2023.9.5
# or
conda install -c conda-forge rdkit
```

### Issue 2: "FileNotFoundError: outputs/models/best_model.pkl"

**Solution**: You haven't trained the model yet
```bash
python src/models.py
```

### Issue 3: API returns "Model not loaded"

**Solution**: Restart the API server
```bash
# Stop the server (Ctrl+C)
# Verify model exists
dir outputs\models\best_model.pkl
# Restart
python -m uvicorn app.backend.main:app --reload --port 8000
```

### Issue 4: Frontend can't connect to API

**Solution**: Check CORS and API URL
```bash
# In app/frontend/.env
echo REACT_APP_API_URL=http://localhost:8000 > .env

# Restart frontend
npm run dev
```

### Issue 5: "Port 8000 already in use"

**Solution**: Kill the process or use a different port
```powershell
# Find process
netstat -ano | findstr :8000

# Kill it (replace PID)
taskkill /PID <PID> /F

# Or use different port
python -m uvicorn app.backend.main:app --reload --port 8001
```

### Issue 6: Tests fail with "Invalid SMILES"

**Solution**: The molecular features module needs SMILES in the dataset
```bash
# Ensure dataset has SMILES column
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_final.csv'); print(df.columns)"

# If SMILES missing, regenerate descriptors
python src/molecular_features.py
```

---

## 📊 Expected Results Summary

After complete setup:

| Component | Status | Location |
|-----------|--------|----------|
| Model trained | ✅ | `outputs/models/best_model.pkl` |
| Preprocessor | ✅ | `outputs/models/preprocessor.pkl` |
| Visualizations | ✅ | `outputs/figures/` (20+ files) |
| Analysis results | ✅ | `outputs/analysis/` (JSON files) |
| Reports | ✅ | `outputs/*.md` (3 reports) |
| API running | ✅ | http://localhost:8000 |
| Frontend running | ✅ | http://localhost:5173 |
| Tests passing | ✅ | 90%+ pass rate |

---

## 🚀 Optional: Advanced Setup

### Enable Production Mode

```bash
# Build frontend for production
cd app/frontend
npm run build

# Serve with backend
cd ../..
python -m uvicorn app.backend.main:app --host 0.0.0.0 --port 8000
```

### Run All Analyses in One Command

Create `run_all.sh` (Git Bash) or `run_all.bat` (Windows):

```batch
@echo off
echo Running complete analysis pipeline...
python src/toxicophores.py
python src/recommendations.py
python src/source_comparison.py
python run_comprehensive_analysis.py
echo Complete! Check outputs/ directory
```

### Set Up Git Hooks (for development)

```bash
# Pre-commit hook for code quality
pip install pre-commit
pre-commit install
```

---

## 📚 Next Steps

1. **Explore the API**: Visit http://localhost:8000/docs for interactive API documentation
2. **Review Reports**: Check `outputs/*.md` for detailed analysis reports
3. **Customize Models**: Edit `src/models.py` to try different algorithms
4. **Add Features**: Modify `src/molecular_features.py` to add new descriptors
5. **Deploy**: Follow `DEPLOYMENT.md` for Vercel/Railway deployment

---

## 🆘 Need Help?

- **API Documentation**: http://localhost:8000/docs
- **Frontend README**: `app/frontend/README.md`
- **Model Card**: `MODEL_CARD.md`
- **Reproducibility**: `REPRODUCIBILITY.md`

---

## ✅ Success Criteria

You've successfully set up the system when:

1. ✅ All dependencies installed without errors
2. ✅ Model trained with >80% accuracy
3. ✅ API returns predictions correctly
4. ✅ Frontend displays results
5. ✅ All major tests pass
6. ✅ 20+ visualizations generated
7. ✅ 3+ analysis reports created

**Congratulations! Your ApisTox system is fully operational!** 🎉

---

*For deployment instructions, see `DEPLOYMENT.md`*
*For development guidelines, see `CONTRIBUTING.md`*




