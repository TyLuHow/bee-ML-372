# ApisTox Critical Files Inventory

**Document Version**: 1.0
**Last Updated**: November 19, 2025
**Purpose**: Comprehensive listing of essential files with purposes and dependencies

---

## Overview

This document catalogs all critical files in the ApisTox project, organized by functional category. Each entry includes:
- **File Path**: Absolute location
- **Size/Lines**: File metrics
- **Purpose**: What the file does
- **Dependencies**: What it requires
- **Used By**: What depends on it
- **Criticality**: Impact if missing (Critical/High/Medium/Low)

---

## 🔴 Critical Files (System Cannot Run Without These)

### Backend API

#### 1. FastAPI Main Application
```
FILE: app/backend/main.py
LINES: 574
SIZE: ~45 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - FastAPI REST API server
  - 10 endpoint definitions
  - Request/response validation (Pydantic)
  - Model and preprocessor loading at startup
  - CORS middleware configuration
  - Prediction history logging

DEPENDENCIES:
  - fastapi, uvicorn, pydantic (pip packages)
  - outputs/models/best_model_xgboost.pkl
  - outputs/preprocessors/preprocessor.pkl
  - src/molecular_features.py (for SMILES predictions)

USED BY:
  - app/frontend/src/services/api.ts (HTTP client)
  - Docker: Dockerfile.backend
  - Deployment: vercel.json, docker-compose.yml

CRITICALITY: 🔴 CRITICAL
  - Without this: API completely non-functional
  - Impact: Frontend cannot make predictions
  - Recovery: No alternative; must fix immediately

ENTRY POINTS:
  - Direct: python app/backend/main.py
  - Uvicorn: uvicorn app.backend.main:app --reload
  - Docker: CMD in Dockerfile.backend

KEY FUNCTIONS:
  - predict(input_data: PredictionInput) → Toxicity prediction
  - predict_from_smiles(input_data: PredictionInputSMILES) → SMILES-based prediction
  - get_model_info() → Model metadata
  - get_feature_importance() → SHAP values
  - get_toxicophore_analysis() → Structural alerts
  - recommend_alternatives(cid: int) → Safer compounds
```

#### 2. Production ML Model
```
FILE: outputs/models/best_model_xgboost.pkl
SIZE: 13.2 MB
FORMAT: Joblib pickle

PURPOSE:
  - Trained XGBoost classifier
  - Binary toxicity prediction (0=non-toxic, 1=toxic)
  - Best performer (83.6% accuracy, 85.8% ROC-AUC)

DEPENDENCIES:
  - xgboost>=2.0.3
  - numpy, pandas, scikit-learn

USED BY:
  - app/backend/main.py (loaded at startup)
  - model.predict(), model.predict_proba()

CRITICALITY: 🔴 CRITICAL
  - Without this: No predictions possible
  - Impact: API /predict endpoints return 500 errors
  - Recovery: Retrain model using src/models.py

GENERATION:
  - Script: src/models.py → ModelTrainer.train_all_models()
  - Training time: ~12 seconds
  - Training data: data/raw/dataset_with_descriptors.csv

HYPERPARAMETERS (tuned):
  - n_estimators: 200
  - max_depth: 8
  - learning_rate: 0.05
  - subsample: 0.8
  - colsample_bytree: 0.8
  - scale_pos_weight: 2.5 (class imbalance)
```

#### 3. Data Preprocessor
```
FILE: outputs/preprocessors/preprocessor.pkl
SIZE: 0.5 MB
FORMAT: Joblib pickle

PURPOSE:
  - DataPreprocessor instance (fitted on training data)
  - StandardScaler for feature normalization
  - One-hot encoder settings (source, toxicity_type)
  - Column names and order enforcement

DEPENDENCIES:
  - scikit-learn>=1.4.0
  - pandas, numpy

USED BY:
  - app/backend/main.py (input transformation)
  - preprocessor.scaler.transform(input_df)

CRITICALITY: 🔴 CRITICAL
  - Without this: Input features not properly scaled
  - Impact: Model predictions completely incorrect
  - Recovery: Retrain and save using src/preprocessing.py

CONTAINS:
  - StandardScaler (mean=0, std=1 fitted on train)
  - Feature names list (24 expected columns)
  - Categorical mappings (source_ECOTOX, source_PPDB, etc.)

GENERATION:
  - Script: src/preprocessing.py → DataPreprocessor.save_preprocessor()
  - Fit on: 724 training samples
  - Excludes: CID, Preferred_name, SMILES, InChI
```

### Frontend Application

#### 4. React Entry Point
```
FILE: app/frontend/src/main.tsx
LINES: 15
SIZE: ~0.5 KB
LANGUAGE: TypeScript + React

PURPOSE:
  - React application bootstrap
  - Mounts App component to #root div
  - Imports global CSS

DEPENDENCIES:
  - react, react-dom (npm packages)
  - app/frontend/src/App.tsx
  - app/frontend/src/index.css

USED BY:
  - app/frontend/index.html (via <script type="module">)
  - Vite build system

CRITICALITY: 🔴 CRITICAL
  - Without this: Frontend does not render
  - Impact: Blank page in browser
  - Recovery: Must fix; no alternative

CODE:
  import React from 'react'
  import ReactDOM from 'react-dom/client'
  import App from './App.tsx'
  import './index.css'

  ReactDOM.createRoot(document.getElementById('root')!).render(
    <React.StrictMode>
      <App />
    </React.StrictMode>,
  )
```

#### 5. React Main Component
```
FILE: app/frontend/src/App.tsx
LINES: 180
SIZE: ~8 KB
LANGUAGE: TypeScript + React

PURPOSE:
  - Root React component
  - State management (result, loading, error)
  - Layout composition (3-column grid)
  - Event handlers (onPrediction, onError)

DEPENDENCIES:
  - react (useState hook)
  - app/frontend/src/components/PredictionForm.tsx
  - app/frontend/src/components/ResultDisplay.tsx
  - app/frontend/src/components/ModelInfo.tsx
  - app/frontend/src/App.css

USED BY:
  - app/frontend/src/main.tsx

CRITICALITY: 🔴 CRITICAL
  - Without this: No UI components render
  - Impact: Application non-functional
  - Recovery: Must fix; central component

STATE MANAGEMENT:
  - result: PredictionResult | null
  - loading: boolean
  - error: string | null

CALLBACKS:
  - handlePrediction(result) → Update result state
  - handleError(error) → Display error message
```

#### 6. API Client Service
```
FILE: app/frontend/src/services/api.ts
LINES: 80
SIZE: ~3 KB
LANGUAGE: TypeScript

PURPOSE:
  - Axios HTTP client for backend communication
  - Type-safe API functions
  - Error handling and logging

DEPENDENCIES:
  - axios (npm package)
  - TypeScript type definitions

USED BY:
  - app/frontend/src/components/PredictionForm.tsx (predictToxicity)
  - app/frontend/src/components/ModelInfo.tsx (getModelInfo)

CRITICALITY: 🔴 CRITICAL
  - Without this: Frontend cannot communicate with backend
  - Impact: All predictions fail
  - Recovery: Must fix; no alternative

CONFIGURATION:
  - Base URL: http://localhost:8000 (development)
  - Environment variable: VITE_API_URL (production)
  - Content-Type: application/json

FUNCTIONS:
  - predictToxicity(input: PredictionInput) → Promise<PredictionResult>
  - getModelInfo() → Promise<ModelInfo>
  - checkHealth() → Promise<{status: string}>
```

### Data

#### 7. Main Dataset
```
FILE: data/raw/dataset_with_descriptors.csv
ROWS: 1,035
COLUMNS: 28
SIZE: ~0.8 MB
FORMAT: CSV (comma-separated)

PURPOSE:
  - ApisTox training dataset
  - 1,035 compounds with 15 molecular descriptors
  - Binary toxicity labels (EPA_binary)
  - Metadata (year, source, chemical type)

DEPENDENCIES:
  - None (raw data)

USED BY:
  - src/preprocessing.py (DataPreprocessor.load_data)
  - src/models.py (training)
  - src/toxicophores.py (structural analysis)
  - src/temporal_analysis.py (time trends)
  - src/chemical_space.py (visualization)

CRITICALITY: 🔴 CRITICAL
  - Without this: Cannot train models
  - Impact: Must use existing trained models only
  - Recovery: Re-download from ECOTOX/PPDB (requires manual curation)

SCHEMA:
  - Identifiers: CID, Preferred_name, SMILES, InChI
  - Metadata: source, year, toxicity_type, chemical_type
  - Flags: insecticide, herbicide, fungicide, other_agrochemical
  - Descriptors: MolecularWeight, LogP, NumHDonors, ... (15 total)
  - Target: EPA_binary (0 or 1)

DATA QUALITY:
  - Missing values: <1%
  - Duplicates: 0 (CID unique)
  - SMILES validity: 100%
  - Temporal range: 1832-2023 (191 years)
```

---

## 🟠 High Priority Files (Required for Full Functionality)

### ML Pipeline Source Code

#### 8. Data Preprocessing Module
```
FILE: src/preprocessing.py
LINES: 593
SIZE: ~45 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - DataPreprocessor class (complete preprocessing pipeline)
  - Feature engineering, scaling, encoding
  - Train/val/test splitting (stratified + scaffold-based)
  - SMOTE imbalance handling
  - Preprocessor persistence

DEPENDENCIES:
  - pandas, numpy, scikit-learn
  - imbalanced-learn (SMOTE)
  - rdkit (scaffold splitting)

USED BY:
  - src/models.py (training pipeline)
  - Direct execution: python src/preprocessing.py

CRITICALITY: 🟠 HIGH
  - Without this: Cannot retrain models or create new preprocessors
  - Impact: Stuck with existing artifacts
  - Recovery: Rewrite preprocessing logic (significant effort)

KEY METHODS:
  - load_data(file_path) → Load CSV
  - prepare_features(exclude_cols) → X, y separation
  - encode_categorical_features(cols) → One-hot encoding
  - split_data(test_size, val_size, stratify) → Stratified split
  - scaffold_split(smiles_col) → Scaffold-based split
  - scale_features(fit=True) → StandardScaler
  - handle_imbalance(method='smote') → Balance classes
  - save_preprocessor(path) → Joblib save
```

#### 9. Model Training Module
```
FILE: src/models.py
LINES: 611
SIZE: ~50 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - ModelTrainer class (train 6 algorithms)
  - Hyperparameter tuning (GridSearchCV, RandomizedSearchCV)
  - Model evaluation (accuracy, F1, ROC-AUC)
  - Model comparison and selection
  - Model persistence

DEPENDENCIES:
  - scikit-learn, xgboost, lightgbm
  - pandas, numpy, joblib

USED BY:
  - Direct execution: python src/models.py
  - Manual retraining workflows

CRITICALITY: 🟠 HIGH
  - Without this: Cannot retrain or improve models
  - Impact: Stuck with current model performance
  - Recovery: Rewrite training logic (moderate effort)

SUPPORTED MODELS (6):
  1. Logistic Regression (baseline)
  2. Random Forest (ensemble)
  3. XGBoost (gradient boosting) ⭐ Best
  4. LightGBM (gradient boosting)
  5. SVM (support vector machine)
  6. MLP (neural network)

KEY METHODS:
  - get_model(model_name) → Instantiate model
  - train_model(X, y, method='grid') → Hyperparameter tuning
  - evaluate_model(X, y) → Metrics calculation
  - train_all_models() → Train and compare all 6
  - save_model(path) → Joblib save
  - load_model(path) → Joblib load
```

#### 10. Molecular Featurization Module
```
FILE: src/molecular_features.py
LINES: 271
SIZE: ~20 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - MolecularFeaturizer class
  - SMILES → 15 RDKit molecular descriptors
  - Batch processing support
  - Error handling for invalid SMILES

DEPENDENCIES:
  - rdkit>=2023.9.5
  - pandas, numpy

USED BY:
  - app/backend/main.py (SMILES predictions)
  - /predict/smiles endpoint

CRITICALITY: 🟠 HIGH
  - Without this: Cannot process SMILES inputs
  - Impact: /predict/smiles endpoint fails
  - Recovery: Manually calculate descriptors or use alternative library

15 DESCRIPTORS:
  • MolecularWeight, LogP, NumHDonors, NumHAcceptors
  • NumRotatableBonds, NumAromaticRings, TPSA
  • NumHeteroatoms, NumRings, NumSaturatedRings
  • NumAliphaticRings, FractionCSP3, MolarRefractivity
  • BertzCT, HeavyAtomCount

KEY METHODS:
  - smiles_to_mol(smiles) → RDKit Mol object
  - calculate_descriptors(mol) → Dict[str, float]
  - smiles_to_descriptors(smiles) → One-step conversion
  - batch_smiles_to_dataframe(smiles_list) → DataFrame
```

### Frontend Components

#### 11. Prediction Form Component
```
FILE: app/frontend/src/components/PredictionForm.tsx
LINES: 221
SIZE: ~12 KB
LANGUAGE: TypeScript + React

PURPOSE:
  - Main input form (28 fields)
  - Form state management (useState)
  - Submit handler (POST /predict)
  - Input validation (min/max ranges)
  - Example data pre-filled

DEPENDENCIES:
  - react (useState hook)
  - app/frontend/src/services/api.ts (predictToxicity)

USED BY:
  - app/frontend/src/App.tsx

CRITICALITY: 🟠 HIGH
  - Without this: Users cannot input data
  - Impact: Frontend unusable for predictions
  - Recovery: Must rebuild form (moderate effort)

FORM SECTIONS:
  1. Compound Information (source, year, toxicity_type)
  2. Chemical Type (insecticide, herbicide, fungicide, other)
  3. Molecular Descriptors (15 fields in grid layout)

STATE: 28 fields (all form inputs)
VALIDATION: HTML5 constraints (min, max, required, step)
```

#### 12. Result Display Component
```
FILE: app/frontend/src/components/ResultDisplay.tsx
LINES: 150
SIZE: ~8 KB
LANGUAGE: TypeScript + React

PURPOSE:
  - Prediction results visualization
  - Color-coded badges (Green=Non-Toxic, Red=Toxic)
  - Confidence meter (circular progress)
  - Probability breakdown (bar chart)
  - Loading states and error messages

DEPENDENCIES:
  - react
  - recharts (optional, for charts)

USED BY:
  - app/frontend/src/App.tsx

CRITICALITY: 🟠 HIGH
  - Without this: Results not displayed
  - Impact: Users don't see predictions
  - Recovery: Create alternative display (moderate effort)

DISPLAYS:
  - prediction_label ("Toxic" or "Non-Toxic")
  - confidence (0.5-1.0 as percentage)
  - probabilities.toxic, probabilities.non_toxic
  - timestamp (ISO 8601)
```

---

## 🟡 Medium Priority Files (Important for Production)

### Configuration

#### 13. Docker Compose Orchestration
```
FILE: docker-compose.yml
LINES: 45
SIZE: ~2 KB
FORMAT: YAML

PURPOSE:
  - Multi-container orchestration
  - Backend + frontend services
  - Network configuration
  - Volume mounts for data persistence

DEPENDENCIES:
  - Dockerfile.backend
  - Dockerfile.frontend
  - docker, docker-compose (system)

USED BY:
  - Deployment: docker-compose up

CRITICALITY: 🟡 MEDIUM
  - Without this: Manual container management required
  - Impact: Harder to deploy full stack
  - Recovery: Use individual docker run commands

SERVICES:
  - backend: Python 3.9, port 8000
  - frontend: Nginx, port 80
  - Volumes: ./outputs (model persistence)
```

#### 14. Backend Dockerfile
```
FILE: Dockerfile.backend
LINES: 25
SIZE: ~1 KB
FORMAT: Dockerfile

PURPOSE:
  - Backend containerization
  - Python 3.9-slim base image
  - Install dependencies (requirements-production.txt)
  - Copy source code and models
  - Expose port 8000, run uvicorn

DEPENDENCIES:
  - requirements-production.txt
  - app/backend/
  - outputs/models/
  - outputs/preprocessors/

USED BY:
  - docker-compose.yml (backend service)
  - Docker Hub / container registry

CRITICALITY: 🟡 MEDIUM
  - Without this: Cannot containerize backend
  - Impact: Must run locally or create alternative Dockerfile
  - Recovery: Recreate Dockerfile (easy)

BUILD:
  docker build -f Dockerfile.backend -t apistox-backend .
RUN:
  docker run -p 8000:8000 apistox-backend
```

#### 15. Frontend Dockerfile
```
FILE: Dockerfile.frontend
LINES: 30
SIZE: ~1.5 KB
FORMAT: Dockerfile

PURPOSE:
  - Frontend containerization (multi-stage build)
  - Stage 1: Node 18, build Vite app
  - Stage 2: Nginx alpine, serve static files

DEPENDENCIES:
  - app/frontend/package.json
  - app/frontend/src/
  - Nginx configuration

USED BY:
  - docker-compose.yml (frontend service)

CRITICALITY: 🟡 MEDIUM
  - Without this: Cannot containerize frontend
  - Impact: Must serve frontend differently
  - Recovery: Alternative web server (Nginx, Apache)

STAGES:
  1. Build: npm install && npm run build
  2. Serve: Copy build/ to /usr/share/nginx/html
```

#### 16. Python Dependencies
```
FILE: requirements.txt
LINES: 35
SIZE: ~1 KB
FORMAT: Text (pip requirements)

PURPOSE:
  - All Python dependencies with pinned versions
  - Development + production packages

DEPENDENCIES:
  - None (defines dependencies for project)

USED BY:
  - pip install -r requirements.txt
  - Dockerfile.backend
  - CI/CD pipelines

CRITICALITY: 🟡 MEDIUM
  - Without this: Manual package installation required
  - Impact: Dependency version mismatches
  - Recovery: Recreate from pip freeze

KEY PACKAGES:
  - fastapi==0.104.1, uvicorn==0.24.0
  - scikit-learn==1.4.0, xgboost==2.0.3
  - rdkit==2023.9.5
  - pandas==2.1.4, numpy==1.26.2
  - shap==0.43.0, lime==0.2.0.1
```

#### 17. Frontend Dependencies
```
FILE: app/frontend/package.json
LINES: 40
SIZE: ~2 KB
FORMAT: JSON

PURPOSE:
  - Node.js dependencies and scripts
  - React, TypeScript, Vite, TailwindCSS

DEPENDENCIES:
  - None (defines dependencies for frontend)

USED BY:
  - npm install
  - npm run dev, npm run build
  - Dockerfile.frontend

CRITICALITY: 🟡 MEDIUM
  - Without this: Cannot install frontend packages
  - Impact: Build fails
  - Recovery: Recreate package.json (moderate effort)

SCRIPTS:
  - dev: vite (development server)
  - build: tsc && vite build (production build)
  - preview: vite preview (preview build)

KEY PACKAGES:
  - react@18.2.0, react-dom@18.2.0
  - typescript@5.2.2
  - vite@5.0.8
  - tailwindcss@3.3.6
  - axios@1.6.0
```

### Analysis Modules

#### 18. Toxicophore Analysis
```
FILE: src/toxicophores.py
LINES: 422
SIZE: ~35 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - ToxicophoreAnalyzer class
  - 20 SMARTS structural alert patterns
  - Statistical testing (Chi-square, enrichment)
  - Visualization (bar charts, scatter plots)

DEPENDENCIES:
  - rdkit (SMARTS matching)
  - scipy (statistical tests)
  - matplotlib, seaborn (plotting)

USED BY:
  - app/backend/main.py (/analysis/toxicophores endpoints)
  - Direct execution: python src/toxicophores.py

CRITICALITY: 🟡 MEDIUM
  - Without this: Toxicophore analysis unavailable
  - Impact: API endpoints return 404
  - Recovery: Use pre-generated results or disable feature

OUTPUTS:
  - outputs/analysis/toxicophore_results.json
  - outputs/figures/toxicophore_enrichment.png
```

#### 19. Compound Recommendations
```
FILE: src/recommendations.py
LINES: 336
SIZE: ~28 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - CompoundRecommender class
  - KNN-based safer alternatives
  - Similarity scoring (Euclidean distance)
  - Batch recommendation generation

DEPENDENCIES:
  - scikit-learn (KNN)
  - pandas, numpy

USED BY:
  - app/backend/main.py (/recommend/alternatives endpoint)
  - Direct execution: python src/recommendations.py

CRITICALITY: 🟡 MEDIUM
  - Without this: Recommendation feature unavailable
  - Impact: API endpoint returns 404
  - Recovery: Use pre-generated alternatives.csv

OUTPUTS:
  - outputs/analysis/alternatives.csv
```

---

## 🟢 Low Priority Files (Optional, Enhances Experience)

### Interpretability

#### 20. SHAP/LIME Interpretability
```
FILE: src/interpretability.py
LINES: 388
SIZE: ~32 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - ModelInterpreter class
  - SHAP TreeExplainer (feature importance)
  - LIME TabularExplainer (local explanations)
  - Visualization (summary plots, waterfall charts)

DEPENDENCIES:
  - shap>=0.43.0 (LARGE package, 50+ MB)
  - lime>=0.2.0.1
  - matplotlib

USED BY:
  - app/backend/main.py (/feature/importance endpoint)
  - Direct execution: python src/interpretability.py

CRITICALITY: 🟢 LOW
  - Without this: Feature importance unavailable
  - Impact: /feature/importance returns empty or error
  - Recovery: Use pre-calculated importance CSV

OUTPUTS:
  - outputs/metrics/feature_importance_shap.csv
  - outputs/figures/shap_summary.png
  - outputs/figures/shap_importance.png
  - outputs/figures/waterfall_*.png

NOTE: SHAP excluded from serverless deployments (size)
```

### Temporal and Spatial Analysis

#### 21. Temporal Analysis
```
FILE: src/temporal_analysis.py
LINES: 354
SIZE: ~30 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - TemporalAnalyzer class
  - Mann-Kendall trend test
  - Decade-by-decade comparison
  - Time series visualization

DEPENDENCIES:
  - scipy (Mann-Kendall test)
  - matplotlib, seaborn

USED BY:
  - Direct execution: python src/temporal_analysis.py

CRITICALITY: 🟢 LOW
  - Without this: Temporal trends unavailable
  - Impact: No API endpoint affected (analysis only)
  - Recovery: Use pre-generated temporal_trends.json

OUTPUTS:
  - outputs/analysis/temporal_trends.json
  - outputs/figures/temporal_trends.png
  - outputs/figures/decade_comparison.png
```

#### 22. Chemical Space Visualization
```
FILE: src/chemical_space.py
LINES: 378
SIZE: ~32 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - ChemicalSpaceVisualizer class
  - PCA (2D/3D projections)
  - t-SNE embedding
  - Interactive Plotly visualizations

DEPENDENCIES:
  - scikit-learn (PCA, t-SNE)
  - plotly (interactive plots)
  - matplotlib

USED BY:
  - Direct execution: python src/chemical_space.py

CRITICALITY: 🟢 LOW
  - Without this: Chemical space viz unavailable
  - Impact: No API endpoint affected (analysis only)
  - Recovery: Use pre-generated PNGs

OUTPUTS:
  - outputs/analysis/chemical_space_results.json
  - outputs/figures/chemical_space_pca.png
  - outputs/figures/chemical_space_tsne.png
```

#### 23. Source Comparison
```
FILE: src/source_comparison.py
LINES: 304
SIZE: ~25 KB
LANGUAGE: Python 3.9+

PURPOSE:
  - SourceComparator class
  - ECOTOX vs PPDB comparison
  - Statistical testing (t-test, Chi-square)
  - Distribution visualization

DEPENDENCIES:
  - scipy (statistical tests)
  - matplotlib, seaborn

USED BY:
  - Direct execution: python src/source_comparison.py

CRITICALITY: 🟢 LOW
  - Without this: Source comparison unavailable
  - Impact: No API endpoint affected (research only)
  - Recovery: Use existing visualizations

OUTPUTS:
  - outputs/figures/source_distribution_*.png
```

---

## Testing Files

#### 24. API Tests
```
FILE: tests/test_api.py
LINES: 250
SIZE: ~18 KB
LANGUAGE: Python (pytest)

PURPOSE:
  - API endpoint unit tests
  - 15 test cases (all endpoints)
  - Request/response validation
  - Error handling tests

DEPENDENCIES:
  - pytest, fastapi.testclient

USED BY:
  - pytest tests/test_api.py
  - CI/CD pipelines

CRITICALITY: 🟢 LOW
  - Without this: No automated API testing
  - Impact: Manual testing required
  - Recovery: Recreate tests (moderate effort)

COVERAGE:
  - Estimated 40% code coverage
  - All endpoints tested
  - Missing: edge cases, performance tests
```

#### 25. Model Tests
```
FILE: tests/test_models.py
LINES: 200
SIZE: ~15 KB
LANGUAGE: Python (pytest)

PURPOSE:
  - Model training unit tests
  - 10 test cases
  - Hyperparameter tuning, metrics, persistence

DEPENDENCIES:
  - pytest, scikit-learn

USED BY:
  - pytest tests/test_models.py

CRITICALITY: 🟢 LOW
  - Without this: No automated training tests
  - Impact: Manual validation required
  - Recovery: Recreate tests
```

#### 26. Preprocessing Tests
```
FILE: tests/test_preprocessing.py
LINES: 180
SIZE: ~14 KB
LANGUAGE: Python (pytest)

PURPOSE:
  - Preprocessing pipeline tests
  - 12 test cases
  - Data loading, scaling, splitting, SMOTE

DEPENDENCIES:
  - pytest, pandas, scikit-learn

USED BY:
  - pytest tests/test_preprocessing.py

CRITICALITY: 🟢 LOW
  - Without this: No automated preprocessing tests
  - Impact: Manual validation required
  - Recovery: Recreate tests
```

---

## Documentation Files

#### 27. Main README
```
FILE: README.md
LINES: 400
SIZE: ~35 KB
FORMAT: Markdown

PURPOSE:
  - Project overview and documentation
  - Quick start guide
  - Installation instructions
  - Usage examples
  - API documentation overview
  - Deployment instructions
  - Citation information

DEPENDENCIES:
  - None

USED BY:
  - GitHub repository (rendered on main page)
  - New users (onboarding)

CRITICALITY: 🟢 LOW
  - Without this: Harder onboarding for new users
  - Impact: Documentation missing
  - Recovery: Recreate from other docs

SECTIONS:
  - Project Overview
  - Features
  - Installation
  - Quick Start
  - Usage (API + Frontend)
  - Model Performance
  - Deployment
  - Citation
```

#### 28. API Documentation
```
FILE: docs/API_DOCS.md
LINES: 500
SIZE: ~42 KB
FORMAT: Markdown

PURPOSE:
  - Comprehensive API reference
  - All 10 endpoints documented
  - Request/response schemas
  - Example curl commands
  - Error codes

DEPENDENCIES:
  - None

USED BY:
  - API users (reference)
  - Integration developers

CRITICALITY: 🟢 LOW
  - Without this: Auto-generated Swagger docs available
  - Impact: Manual documentation missing
  - Recovery: Use /docs (FastAPI auto-docs)
```

#### 29. Model Card
```
FILE: docs/MODEL_CARD.md
LINES: 300
SIZE: ~25 KB
FORMAT: Markdown

PURPOSE:
  - ML model documentation (best practice)
  - Model details, intended use, limitations
  - Ethical considerations
  - Performance metrics

DEPENDENCIES:
  - None

USED BY:
  - Regulatory submissions
  - Research papers
  - Transparency reporting

CRITICALITY: 🟢 LOW
  - Without this: Less transparency
  - Impact: Documentation gap
  - Recovery: Recreate from training results
```

---

## Generated Artifacts (Auto-Created)

### Model Artifacts
```
outputs/models/
  - best_model_xgboost.pkl        [13.2 MB] CRITICAL
  - best_model_random_forest.pkl  [8.7 MB]  HIGH
  - best_model_lightgbm.pkl       [2.1 MB]  MEDIUM
  - best_model_svm.pkl            [5.3 MB]  MEDIUM
  - best_model_mlp.pkl            [1.8 MB]  MEDIUM
  - best_model_logistic.pkl       [0.3 MB]  LOW

Generated by: src/models.py → ModelTrainer.train_all_models()
Regeneration time: ~5 minutes (all models)
Dependencies: data/raw/dataset_with_descriptors.csv
```

### Preprocessor Artifacts
```
outputs/preprocessors/
  - preprocessor.pkl              [0.5 MB]  CRITICAL

Generated by: src/preprocessing.py → DataPreprocessor.save_preprocessor()
Regeneration time: ~30 seconds
Dependencies: data/raw/dataset_with_descriptors.csv
```

### Metrics
```
outputs/metrics/
  - training_results.json         [8 KB]    HIGH
  - feature_importance_shap.csv   [2 KB]    MEDIUM
  - confusion_matrix.csv          [1 KB]    LOW
  - classification_report.txt     [1 KB]    LOW

Generated by: src/models.py, src/interpretability.py
Regeneration time: ~2 minutes
```

### Visualizations
```
outputs/figures/
  - shap_summary.png              [1.2 MB]  MEDIUM
  - shap_importance.png           [0.8 MB]  MEDIUM
  - toxicophore_enrichment.png    [0.9 MB]  MEDIUM
  - chemical_space_pca.png        [1.5 MB]  LOW
  - temporal_trends.png           [1.1 MB]  LOW
  - confusion_matrix_plot.png     [0.6 MB]  LOW
  - roc_curve.png                 [0.7 MB]  LOW
  - (13+ additional plots)

Generated by: src/*.py analysis scripts
Regeneration time: ~5 minutes (all plots)
```

### Analysis Results
```
outputs/analysis/
  - toxicophore_results.json      [12 KB]   MEDIUM
  - alternatives.csv              [45 KB]   MEDIUM
  - temporal_trends.json          [8 KB]    LOW
  - chemical_space_results.json   [15 KB]   LOW

Generated by: src/toxicophores.py, src/recommendations.py, etc.
Regeneration time: ~3 minutes
```

---

## File Dependency Graph

```
┌─────────────────────────────────────────────────────────────┐
│ DATA LAYER                                                  │
│ dataset_with_descriptors.csv                                │
└────────────────────────┬────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ PREPROCESSING LAYER                                         │
│ src/preprocessing.py → preprocessor.pkl                     │
└────────────────────────┬────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ TRAINING LAYER                                              │
│ src/models.py → best_model_xgboost.pkl                      │
└────────────────────────┬────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ API LAYER                                                   │
│ app/backend/main.py (loads model + preprocessor)            │
└────────────────────────┬────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────┐
│ FRONTEND LAYER                                              │
│ app/frontend/src/services/api.ts → App.tsx                  │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ PARALLEL: ANALYSIS LAYER (Optional)                         │
│ src/molecular_features.py → SMILES predictions              │
│ src/toxicophores.py → toxicophore_results.json              │
│ src/recommendations.py → alternatives.csv                   │
│ src/interpretability.py → feature_importance_shap.csv       │
│ src/temporal_analysis.py → temporal_trends.json             │
│ src/chemical_space.py → chemical_space_results.json         │
└─────────────────────────────────────────────────────────────┘
```

---

## Criticality Matrix

| Criticality | Count | Total Size | Impact if Missing |
|-------------|-------|------------|-------------------|
| 🔴 **CRITICAL** | 7 | ~60 MB | System non-functional |
| 🟠 **HIGH** | 12 | ~5 MB | Major features broken |
| 🟡 **MEDIUM** | 15 | ~10 MB | Production deployment affected |
| 🟢 **LOW** | 20+ | ~100 MB | Enhanced features unavailable |

---

## Recovery Procedures

### If Critical Files Missing

**1. Model Missing (best_model_xgboost.pkl)**:
```bash
# Regenerate model (5 minutes)
python src/preprocessing.py  # Creates preprocessor.pkl
python src/models.py          # Creates all models
```

**2. Preprocessor Missing (preprocessor.pkl)**:
```bash
# Regenerate preprocessor (30 seconds)
python src/preprocessing.py
```

**3. Dataset Missing (dataset_with_descriptors.csv)**:
- **Cannot auto-regenerate** (requires manual data curation)
- Recovery: Download from ECOTOX/PPDB and reprocess (days of work)
- Backup: Keep multiple copies in separate locations

**4. API Server Missing (app/backend/main.py)**:
- **Cannot auto-regenerate** (requires rewriting)
- Recovery: Restore from git history
- Backup: Version control essential

**5. Frontend Components Missing**:
- **Cannot auto-regenerate** (requires rewriting UI)
- Recovery: Restore from git history
- Backup: Version control essential

### If High Priority Files Missing

**Preprocessing/Model Scripts**:
- Impact: Cannot retrain models
- Workaround: Use existing trained models
- Recovery: Rewrite from scratch (days of work)

**Frontend Components**:
- Impact: Degraded UI
- Workaround: Use alternative input methods (curl, Postman)
- Recovery: Recreate components (hours of work)

### If Medium/Low Priority Files Missing

**Analysis Modules**:
- Impact: Some API endpoints return errors
- Workaround: Use pre-generated results
- Recovery: Rewrite if needed (hours to days)

**Documentation**:
- Impact: Onboarding harder
- Workaround: Use auto-generated docs (/docs)
- Recovery: Recreate from code (hours)

---

## Backup Recommendations

### Critical (Daily Backups)
- `data/raw/dataset_with_descriptors.csv`
- `outputs/models/best_model_xgboost.pkl`
- `outputs/preprocessors/preprocessor.pkl`
- `app/backend/main.py`
- `app/frontend/src/` (entire directory)

### High Priority (Weekly Backups)
- `src/` (all Python source code)
- `requirements.txt`
- `docker-compose.yml`, Dockerfiles

### Medium/Low Priority (Monthly Backups or Git Only)
- Generated artifacts (can be regenerated)
- Documentation (tracked in git)
- Tests (tracked in git)

---

## File Size Breakdown

### By Category
```
Models (6 files):           ~31 MB  (68%)
Figures (20+ files):        ~12 MB  (26%)
Data (1 file):              ~0.8 MB (2%)
Code (25 files):            ~0.5 MB (1%)
Config (15 files):          ~0.2 MB (<1%)
Docs (8 files):             ~0.3 MB (1%)
Tests (3 files):            ~0.1 MB (<1%)
Other:                      ~0.5 MB (1%)
────────────────────────────────────────
TOTAL (excl. node_modules): ~45 MB
```

### Top 10 Largest Files
```
1.  best_model_xgboost.pkl        13.2 MB
2.  best_model_random_forest.pkl   8.7 MB
3.  best_model_svm.pkl             5.3 MB
4.  best_model_lightgbm.pkl        2.1 MB
5.  best_model_mlp.pkl             1.8 MB
6.  chemical_space_pca.png         1.5 MB
7.  chemical_space_tsne.png        1.3 MB
8.  shap_summary.png               1.2 MB
9.  temporal_trends.png            1.1 MB
10. correlation_matrix.png         2.0 MB
```

---

**Document Maintained By**: ApisTox Development Team
**Last Review**: November 19, 2025
**Next Review**: December 2025
