# Honey Bee Toxicity Prediction System
## IME 372 Course Project - Predictive Analytics

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://python.org)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.104%2B-green)](https://fastapi.tiangolo.com)
[![XGBoost](https://img.shields.io/badge/XGBoost-2.0%2B-orange)](https://xgboost.readthedocs.io)
[![License](https://img.shields.io/badge/License-CC--BY--NC--4.0-yellow)](https://creativecommons.org/licenses/by-nc/4.0/)

**A comprehensive machine learning system for predicting pesticide toxicity to honey bees using molecular descriptors and agrochemical properties.**

---

## 🎯 Project Overview

This project implements an end-to-end machine learning pipeline to predict whether pesticides are toxic to honey bees, addressing a critical agricultural and environmental challenge. The system achieves **83.6% accuracy** and **85.8% ROC-AUC** on the test set.

### Key Features

#### Machine Learning & Prediction
- ✅ **Multiple ML Models**: Logistic Regression, Random Forest, XGBoost, LightGBM
- ✅ **Model Interpretability**: SHAP and LIME explanations
- ✅ **SMILES Support**: Direct prediction from molecular structures
- ✅ **Feature Engineering**: 15 molecular descriptors from SMILES
- ✅ **Class Imbalance Handling**: SMOTE resampling

#### Data Explorer (NEW)
- ✅ **Dataset Overview**: Comprehensive statistics and distributions
- ✅ **Molecular Diversity**: Descriptor distributions with toxicity overlay
- ✅ **Toxicity Analysis**: By chemical class with statistical tests
- ✅ **Temporal Trends**: Historical toxicity patterns with Mann-Kendall test
- ✅ **Chemical Space**: PCA and t-SNE visualizations
- ✅ **Toxicophores**: Structural alert enrichment analysis
- ✅ **Feature Correlations**: Interactive correlation matrix
- ✅ **Property Relationships**: 2D scatter plots

#### Web Application (NEW)
- ✅ **Modern Dashboard**: Overview with quick actions
- ✅ **Interactive Explorer**: 7 tabs for comprehensive data analysis
- ✅ **Enhanced Prediction**: SMILES input with real-time validation
- ✅ **Visual Results**: Confidence gauges and toxicity indicators
- ✅ **Responsive Design**: Optimized for desktop and mobile
- ✅ **Performance**: Code splitting and caching for fast load times

#### Backend Infrastructure
- ✅ **REST API**: FastAPI with automatic OpenAPI documentation
- ✅ **Caching System**: In-memory caching for expensive computations
- ✅ **Optimized**: Reduced code duplication with utility modules
- ✅ **Production Ready**: Error handling and validation

---

## 📊 Project Results

### Model Performance (Test Set)

| Metric | Score |
|--------|-------|
| **Accuracy** | 83.57% |
| **F1 Score** | 70.18% |
| **ROC-AUC** | 85.83% |
| **Precision (Toxic)** | 72.73% |
| **Recall (Toxic)** | 67.80% |

### Top Predictive Features (SHAP Analysis)

1. **Insecticide** (1.366) - Chemical type flag
2. **Herbicide** (1.054) - Chemical type flag
3. **Fungicide** (0.740) - Chemical type flag
4. **Year** (0.641) - Publication year
5. **LogP** (0.474) - Lipophilicity

---

## 🗂️ Project Structure

```
apis_tox_dataset/
├── app/
│   └── backend/
│       └── main.py                 # FastAPI REST API
├── data/
│   └── raw/
│       └── dataset_with_descriptors.csv
├── notebooks/
│   └── 01_exploratory_analysis.ipynb  # EDA notebook
├── src/
│   ├── preprocessing.py            # Data preprocessing pipeline
│   ├── models.py                   # Model training and evaluation
│   └── interpretability.py         # SHAP/LIME analysis
├── outputs/
│   ├── models/                     # Trained model files
│   ├── preprocessors/              # Saved preprocessors
│   ├── figures/                    # Generated visualizations
│   └── metrics/                    # Performance metrics
├── tests/
│   └── test_*.py                   # Unit tests
├── docs/
│   ├── project_proposal.md         # Project proposal
│   └── presentation/               # Presentation materials
├── requirements.txt                # Python dependencies
└── README.md                       # This file
```

---

## 🚀 Getting Started

### Prerequisites

- **Backend**: Python 3.8+
- **Frontend**: Node.js 16+ and npm
- **RAM**: 4GB minimum
- **Storage**: 500MB for dependencies

### Installation

#### 1. Clone the Repository
```bash
git clone <repository-url>
cd bee-ML-372
```

#### 2. Backend Setup
```bash
# Install Python dependencies
pip install -r requirements.txt
```

#### 3. Frontend Setup
```bash
# Navigate to frontend directory
cd app/frontend

# Install Node dependencies
npm install

# Return to root
cd ../..
```

### Quick Start

#### Option 1: Run Full Application

**Terminal 1 - Backend**:
```bash
cd /home/user/bee-ML-372
python -m uvicorn app.backend.main:app --host 0.0.0.0 --port 8000 --reload
```

**Terminal 2 - Frontend**:
```bash
cd app/frontend
npm run dev
```

The application will be available at:
- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8000
- **API Docs**: http://localhost:8000/docs

#### Option 2: API Only

```bash
python app/backend/main.py
```

Visit http://localhost:8000/docs for interactive API documentation.

### First Steps

1. **Explore the Dataset**: Navigate to the Explorer tab to analyze the dataset
2. **Make Predictions**: Go to Predict tab and enter a SMILES string or molecular descriptors
3. **View Model Info**: Check the Model tab for performance metrics
4. **API Integration**: Use the /docs endpoint to test API calls

---

## 📈 Dataset Information

**Source**: ApisTox Dataset  
**Size**: 1,035 pesticide compounds  
**Target**: Binary classification (0=non-toxic, 1=toxic)  
**Features**: 24 total (after preprocessing)
- 15 molecular descriptors (from SMILES)
- 7 agrochemical flags
- 2 temporal/source features

**Class Distribution**:
- Non-toxic: 739 (71.4%)
- Toxic: 296 (28.6%)
- Imbalance ratio: 2.50:1

---

## 🔬 Methodology

### 1. Data Preprocessing
- ✅ No missing values
- ✅ Molecular descriptor extraction using RDKit
- ✅ One-hot encoding for categorical features
- ✅ StandardScaler for numerical features
- ✅ Stratified train/val/test split (70/10/20)
- ✅ SMOTE resampling for class imbalance

### 2. Model Development
- **Baseline**: Logistic Regression
- **Ensemble**: Random Forest, XGBoost, LightGBM
- **SVM**: Support Vector Machine
- **Neural Network**: Multi-Layer Perceptron

**Best Model**: XGBoost Classifier with hyperparameter tuning

### 3. Model Interpretability
- **SHAP**: Global and local feature importance
- **LIME**: Individual prediction explanations
- **Feature Analysis**: Correlation and dependency plots

### 4. Production Deployment
- **FastAPI**: RESTful API with automatic documentation
- **Model Serving**: Joblib persistence
- **Prediction History**: SQLite/JSON storage
- **CORS**: Frontend integration ready

---

## 🔌 API Endpoints

### Core Endpoints

#### Health Check
```bash
GET /health
```

#### Make Prediction (Descriptors)
```bash
POST /predict
Content-Type: application/json

{
  "source": "PPDB",
  "year": 2020,
  "toxicity_type": "Contact",
  "insecticide": 1,
  "MolecularWeight": 350.0,
  "LogP": 3.5,
  ...
}
```

#### Make Prediction (SMILES)
```bash
POST /predict/smiles
Content-Type: application/json

{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "year": 2024,
  "insecticide": 1,
  "source": "PPDB",
  "toxicity_type": "Contact"
}
```

#### Model Information
```bash
GET /model/info
```

#### Feature Importance
```bash
GET /feature/importance
```

#### Prediction History
```bash
GET /history?limit=10
```

### Data Explorer Endpoints

#### Dataset Overview
```bash
GET /api/explorer/overview
```

#### Molecular Diversity
```bash
GET /api/explorer/molecular-diversity
```

#### Toxicity by Class
```bash
GET /api/explorer/toxicity-by-class
```

#### Temporal Trends
```bash
GET /api/explorer/temporal-trends
```

#### Chemical Space (PCA/t-SNE)
```bash
GET /api/explorer/chemical-space
```

#### Toxicophores
```bash
GET /api/explorer/toxicophores
```

#### Feature Correlations
```bash
GET /api/explorer/correlations
```

#### Property Distributions
```bash
GET /api/explorer/property-distributions
```

See full API documentation at `/docs` endpoint.

---

## 📊 Visualizations

The project generates comprehensive visualizations:

1. **EDA Visualizations**:
   - Target distribution (bar and pie charts)
   - Molecular descriptor distributions
   - Feature correlations heatmap
   - Toxicity comparison boxplots

2. **Model Interpretability**:
   - SHAP summary plots (beeswarm)
   - SHAP feature importance (bar chart)
   - SHAP waterfall plots (individual predictions)
   - LIME explanations

All visualizations saved in `outputs/figures/`

---

## 🧪 Testing

Run unit tests:
```bash
pytest tests/ -v
```

Test API endpoints:
```bash
python test_api.py
```

---

## 📚 Academic Deliverables

### 1. Project Proposal (2-3 pages)
- **Location**: `docs/project_proposal.md`
- **Contents**: Problem statement, methodology, timeline, team roles

### 2. Presentation Materials (12-15 minutes)
- **Location**: `docs/presentation/`
- **Contents**: Slides covering all project phases with visualizations

### 3. Technical Documentation
- **README.md**: Project overview and usage guide
- **Notebooks**: Interactive Jupyter analysis
- **Code Comments**: Comprehensive docstrings

---

## 🌍 Ethical Considerations

### Environmental Impact
- **Pollinator Conservation**: Models inform safe pesticide use
- **Bee Population Health**: Predictions support regulatory decisions
- **Agricultural Sustainability**: Balance crop protection with bee safety

### Model Limitations
- **Data Bias**: Dataset may not represent all pesticide types
- **Probabilistic**: Predictions are not definitive assessments
- **Transparency**: Full interpretability provided via SHAP/LIME
- **Precautionary Principle**: Uncertainty should favor bee safety

### Responsible Use
- Tool aids conservation, not harmful development
- Stakeholders: farmers, beekeepers, regulators, manufacturers
- Data lineage: public chemical data (no privacy concerns)

---

## 🏆 Key Achievements

✅ **Data Analysis**: Comprehensive EDA with 1,035 pesticide compounds  
✅ **Feature Engineering**: 15 molecular descriptors extracted from SMILES  
✅ **Model Training**: 4 ML algorithms trained and compared  
✅ **Performance**: 83.6% accuracy, 85.8% ROC-AUC  
✅ **Interpretability**: SHAP and LIME explanations implemented  
✅ **API Development**: Production-ready FastAPI backend  
✅ **Documentation**: Academic-grade reports and visualizations  
✅ **Reproducibility**: All results documented with random seeds  

---

## 📖 References

1. **ApisTox Dataset**: [Scientific Data (2024)](https://www.nature.com/articles/s41597-024-04232-w)
2. **SHAP**: Lundberg & Lee (2017). "A Unified Approach to Interpreting Model Predictions"
3. **XGBoost**: Chen & Guestrin (2016). "XGBoost: A Scalable Tree Boosting System"
4. **RDKit**: Open-source cheminformatics toolkit
5. **SMOTE**: Chawla et al. (2002). "SMOTE: Synthetic Minority Over-sampling Technique"

---

## 👥 Team Information

**Course**: IME 372 - Predictive Analytics  
**Institution**: [University Name]  
**Semester**: Fall 2025  
**Project Type**: Machine Learning Classification with Interpretability  

---

## 📝 License

This project uses the ApisTox dataset, which is licensed under [CC-BY-NC-4.0](https://creativecommons.org/licenses/by-nc/4.0/).

**Non-Commercial Use Only**: This project is for educational and research purposes.

---

## 🆘 Support

For questions or issues:
1. Check the documentation in `docs/`
2. Review API documentation at `/docs` endpoint
3. Examine example notebooks in `notebooks/`
4. Contact project team

---

## 🎓 Acknowledgments

- **ApisTox Team**: For the comprehensive dataset
- **Scientific Community**: For RDKit, SHAP, and ML libraries
- **Course Instructors**: For project guidance and support

---

## 📈 Completed Enhancements (v2.0.0)

- ✅ **Frontend**: React/TypeScript web application with Tailwind CSS
- ✅ **Data Explorer**: 8 comprehensive endpoints for dataset analysis
- ✅ **Advanced Visualizations**: Interactive charts with Recharts
- ✅ **Code Splitting**: Lazy loading for optimized bundle size
- ✅ **Caching System**: In-memory caching for performance
- ✅ **SMILES Support**: Direct molecular structure input
- ✅ **Responsive Design**: Mobile and desktop optimized
- ✅ **Production Ready**: Error handling, validation, and optimization

## 🚧 Future Enhancements

- [ ] **Database**: PostgreSQL for production storage
- [ ] **Monitoring**: MLflow for experiment tracking
- [ ] **Docker**: Full containerization with docker-compose
- [ ] **CI/CD**: Automated testing and deployment pipeline
- [ ] **Deep Learning**: Graph Neural Networks for molecular prediction
- [ ] **User Authentication**: Secure user accounts and sessions
- [ ] **Batch Predictions**: Upload CSV for bulk processing
- [ ] **Mobile App**: iOS/Android applications

---

**Built with ❤️ for honey bees and sustainable agriculture**
