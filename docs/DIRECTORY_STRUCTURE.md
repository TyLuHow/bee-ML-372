# ApisTox Directory Structure

**Document Version**: 1.0
**Last Updated**: November 19, 2025
**Purpose**: Complete annotated directory tree with file descriptions

---

## Root Directory Overview

```
/home/user/bee-ML-372/
├── Root configuration files
├── app/                    # Application layer (frontend + backend)
├── src/                    # ML pipeline source code
├── data/                   # Datasets and raw data
├── outputs/                # Generated artifacts (models, figures, metrics)
├── tests/                  # Unit and integration tests
├── docs/                   # Documentation files
└── scripts/                # Utility and automation scripts
```

---

## Complete Directory Tree

### 📁 Root Level

```
/home/user/bee-ML-372/
│
├── 📄 README.md                           [400 lines]
│   └── Comprehensive project documentation
│       • Quick start guide
│       • Installation instructions
│       • Usage examples
│       • API documentation overview
│       • Deployment instructions
│       • Citation information
│
├── 📄 START_HERE.md                       [80 lines]
│   └── Quick start guide for new users
│       • 5-minute setup guide
│       • Essential commands
│       • First prediction example
│
├── 📄 PROJECT_SUMMARY.md                  [120 lines]
│   └── High-level project overview
│       • Problem statement
│       • Solution approach
│       • Key features
│       • Results summary
│
├── 📄 requirements.txt                    [35 lines]
│   └── Python dependencies (full)
│       • FastAPI, scikit-learn, XGBoost
│       • RDKit, SHAP, pandas, numpy
│       • Matplotlib, seaborn, plotly
│       • Complete with pinned versions
│
├── 📄 requirements-production.txt         [12 lines]
│   └── Production-only dependencies
│       • FastAPI + uvicorn
│       • Core ML libraries
│       • Excludes development tools
│       • Optimized for serverless
│
├── 📄 requirements-vercel.txt             [10 lines]
│   └── Vercel serverless dependencies
│       • Minimal subset for deployment
│       • Compatible with Python 3.9
│       • Numpy/Pandas version constraints
│
├── 📄 docker-compose.yml                  [45 lines]
│   └── Multi-container orchestration
│       • Backend service (Python 3.9)
│       • Frontend service (Node 18)
│       • Network configuration
│       • Volume mounts
│       • Environment variables
│
├── 📄 Dockerfile.backend                  [25 lines]
│   └── Backend containerization
│       • Base: python:3.9-slim
│       • Install dependencies
│       • Copy source code
│       • Expose port 8000
│       • Run uvicorn server
│
├── 📄 Dockerfile.frontend                 [30 lines]
│   └── Frontend containerization
│       • Base: node:18 (build stage)
│       • Build Vite app
│       • Base: nginx:alpine (serve stage)
│       • Copy built assets
│       • Nginx configuration
│
├── 📄 vercel.json                         [15 lines]
│   └── Vercel serverless deployment config
│       • Route API to Python functions
│       • Environment variables
│       • Build settings
│
├── 📄 .gitignore                          [80 lines]
│   └── Git ignore patterns
│       • Python cache (__pycache__, *.pyc)
│       • Virtual environments (venv/)
│       • Node modules (node_modules/)
│       • Build artifacts (dist/, build/)
│       • IDE files (.vscode/, .idea/)
│
└── 📄 .env.example                        [10 lines]
    └── Environment variable template
        • API_URL configuration
        • PORT settings
        • Debug flags
        • Secret key placeholders
```

---

### 📁 app/ - Application Layer

```
app/
├── backend/                               # FastAPI REST API
│   ├── 📄 main.py                        [574 lines] ⭐ CRITICAL
│   │   └── FastAPI application entry point
│   │       • 10 RESTful endpoints
│   │       • Model/preprocessor loading
│   │       • Pydantic request/response models
│   │       • CORS middleware configuration
│   │       • Error handling and validation
│   │       • Prediction history logging
│   │       • Automatic OpenAPI docs (/docs, /redoc)
│   │       ENDPOINTS:
│   │       • GET  /              → API information
│   │       • GET  /health        → Health check
│   │       • POST /predict       → Toxicity prediction (descriptors)
│   │       • POST /predict/smiles → Toxicity prediction (SMILES)
│   │       • GET  /model/info    → Model metadata
│   │       • GET  /feature/importance → SHAP values
│   │       • GET  /history       → Prediction log
│   │       • GET  /analysis/toxicophores → Toxicophore stats
│   │       • POST /analysis/toxicophores/molecule → Molecule analysis
│   │       • GET  /recommend/alternatives/{cid} → Safer alternatives
│   │
│   ├── 📄 prediction_history.json        [Auto-generated]
│   │   └── Logged predictions (last 100)
│   │       • Timestamp, input features, prediction
│   │       • Used for analytics and debugging
│   │
│   └── 📄 __init__.py                    [Empty]
│       └── Python package marker
│
└── frontend/                              # React Web Application
    ├── 📄 package.json                   [40 lines]
    │   └── Node.js dependencies and scripts
    │       • react@18.2.0, typescript@5.2.2
    │       • vite@5.0.8 (build tool)
    │       • tailwindcss@3.3.6 (styling)
    │       • axios@1.6.0 (HTTP client)
    │       • recharts@2.10.0 (charts)
    │       SCRIPTS:
    │       • npm run dev     → Development server (port 5173)
    │       • npm run build   → Production build
    │       • npm run preview → Preview build
    │
    ├── 📄 package-lock.json              [Auto-generated]
    │   └── Locked dependency versions
    │
    ├── 📄 vite.config.ts                 [15 lines]
    │   └── Vite build configuration
    │       • React plugin
    │       • Port 5173
    │       • Proxy to backend (/api → :8000)
    │
    ├── 📄 tsconfig.json                  [25 lines]
    │   └── TypeScript compiler configuration
    │       • Target: ES2020
    │       • JSX: react-jsx
    │       • Strict mode enabled
    │
    ├── 📄 tsconfig.node.json             [10 lines]
    │   └── TypeScript config for Node.js
    │
    ├── 📄 tailwind.config.js             [20 lines]
    │   └── TailwindCSS configuration
    │       • Custom color palette
    │       • Gradient utilities
    │       • Responsive breakpoints
    │
    ├── 📄 postcss.config.js              [8 lines]
    │   └── PostCSS configuration
    │       • Tailwind plugin
    │       • Autoprefixer
    │
    ├── 📄 index.html                     [20 lines]
    │   └── HTML entry point
    │       • Root div (#root)
    │       • Loads main.tsx
    │       • Meta tags, title
    │
    ├── 📄 .eslintrc.cjs                  [30 lines]
    │   └── ESLint configuration
    │       • TypeScript rules
    │       • React hooks rules
    │       • Code quality enforcement
    │
    ├── public/                           # Static assets
    │   └── 📄 vite.svg                  [Logo]
    │
    └── src/                              # React source code
        ├── 📄 main.tsx                   [15 lines] ⭐ ENTRY POINT
        │   └── React application entry point
        │       • Renders App component
        │       • Mounts to #root
        │       • Imports global CSS
        │
        ├── 📄 App.tsx                    [180 lines] ⭐ MAIN COMPONENT
        │   └── Root React component
        │       • State management (result, loading, error)
        │       • Layout: 3-column grid
        │       • Gradient background
        │       • Callback handling (onPrediction, onError)
        │       STRUCTURE:
        │       • Left: PredictionForm
        │       • Right-Top: ResultDisplay
        │       • Right-Bottom: ModelInfo
        │
        ├── 📄 App.css                    [50 lines]
        │   └── Component-specific styles
        │       • Gradient animations
        │       • Card shadows
        │       • Responsive utilities
        │
        ├── 📄 index.css                  [30 lines]
        │   └── Global styles
        │       • Tailwind imports
        │       • CSS reset
        │       • Base typography
        │
        ├── components/                   # React components
        │   ├── 📄 PredictionForm.tsx    [221 lines] ⭐ MAIN FORM
        │   │   └── Compound input form
        │   │       • 28 form fields (state management)
        │   │       • Input validation
        │   │       • Submit handler (axios POST)
        │   │       • Example data pre-filled
        │   │       • Grouped sections:
        │   │         - Compound Information (source, year, type)
        │   │         - Chemical Type (insecticide, herbicide, etc.)
        │   │         - Molecular Descriptors (15 fields in grid)
        │   │       • Loading state (disabled inputs)
        │   │       • Error display
        │   │
        │   ├── 📄 ResultDisplay.tsx     [150 lines]
        │   │   └── Prediction results visualization
        │   │       • Color-coded badge (Green/Red)
        │   │       • Confidence meter (circular progress)
        │   │       • Probability breakdown (bars)
        │   │       • Timestamp display
        │   │       • Loading spinner
        │   │       • Error messages
        │   │
        │   └── 📄 ModelInfo.tsx         [80 lines]
        │       └── Model metadata display
        │           • Algorithm name (XGBoost)
        │           • Feature count (24)
        │           • Performance metrics (accuracy, ROC-AUC)
        │           • Training date
        │           • Fetches from /model/info API
        │
        └── services/                     # API client layer
            └── 📄 api.ts                [80 lines] ⭐ API CLIENT
                └── Axios HTTP client
                    • Base URL: localhost:8000 (configurable)
                    • Content-Type: application/json
                    TYPE DEFINITIONS:
                    • PredictionInput (40+ fields)
                    • PredictionResult
                    • ModelInfo
                    FUNCTIONS:
                    • predictToxicity(input) → POST /predict
                    • getModelInfo() → GET /model/info
                    • checkHealth() → GET /health
                    • Error handling with try/catch
```

---

### 📁 src/ - ML Pipeline Core

```
src/
├── 📄 preprocessing.py                    [593 lines] ⭐ CRITICAL
│   └── Data preprocessing pipeline
│       CLASS: DataPreprocessor
│       METHODS:
│       • load_data(file_path) → Load CSV
│       • prepare_features(exclude_cols) → X, y separation
│       • encode_categorical_features(cols) → One-hot encoding
│       • split_data(test_size, val_size, stratify) → Train/val/test
│       • scaffold_split(smiles_col) → Scaffold-based splitting
│       • scale_features(fit) → StandardScaler
│       • handle_imbalance(method='smote') → SMOTE/undersampling
│       • save_preprocessor(path) → Joblib persistence
│       • load_preprocessor(path) → Load saved instance
│       PIPELINE:
│       1. Load CSV → DataFrame
│       2. Separate features (X) and target (y)
│       3. One-hot encode categorical (source, toxicity_type)
│       4. Stratified split (70/10/20)
│       5. Fit StandardScaler on train
│       6. Apply SMOTE to train (balance classes)
│       7. Save preprocessor instance
│
├── 📄 models.py                           [611 lines] ⭐ CRITICAL
│   └── Model training and evaluation
│       CLASS: ModelTrainer
│       SUPPORTED MODELS (6):
│       • Logistic Regression (baseline)
│       • Random Forest (ensemble)
│       • XGBoost (gradient boosting) ⭐ Best
│       • LightGBM (gradient boosting)
│       • SVM (support vector machine)
│       • MLP (neural network)
│       METHODS:
│       • get_model(model_name) → Instantiate model
│       • train_model(X, y, method='grid') → Hyperparameter tuning
│       • evaluate_model(X, y) → Metrics calculation
│       • train_all_models() → Train and compare all 6
│       • save_model(path) → Joblib persistence
│       • load_model(path) → Load saved model
│       HYPERPARAMETER TUNING:
│       • GridSearchCV (exhaustive)
│       • RandomizedSearchCV (faster)
│       • 5-fold stratified cross-validation
│       EVALUATION METRICS:
│       • Accuracy, Precision, Recall, F1
│       • ROC-AUC, Confusion Matrix
│       • Classification Report
│       • Training time
│       OUTPUT:
│       • Best model saved to outputs/models/
│       • Metrics saved to outputs/metrics/
│
├── 📄 molecular_features.py               [271 lines] ⭐ CRITICAL
│   └── RDKit molecular descriptor calculation
│       CLASS: MolecularFeaturizer
│       15 DESCRIPTORS CALCULATED:
│       • MolecularWeight → Descriptors.MolWt
│       • LogP → Crippen.MolLogP
│       • NumHDonors → Lipinski.NumHDonors
│       • NumHAcceptors → Lipinski.NumHAcceptors
│       • NumRotatableBonds → Lipinski.NumRotatableBonds
│       • NumAromaticRings → Lipinski.NumAromaticRings
│       • TPSA → rdMolDescriptors.CalcTPSA
│       • NumHeteroatoms → Lipinski.NumHeteroatoms
│       • NumRings → Lipinski.RingCount
│       • NumSaturatedRings → Lipinski.NumSaturatedRings
│       • NumAliphaticRings → Lipinski.NumAliphaticRings
│       • FractionCSP3 → Lipinski.FractionCSP3
│       • MolarRefractivity → Crippen.MolMR
│       • BertzCT → Descriptors.BertzCT
│       • HeavyAtomCount → Lipinski.HeavyAtomCount
│       METHODS:
│       • smiles_to_mol(smiles) → RDKit Mol object
│       • calculate_descriptors(mol) → Dict of 15 descriptors
│       • smiles_to_descriptors(smiles) → One-step conversion
│       • batch_smiles_to_dataframe(smiles_list) → DataFrame
│       ERROR HANDLING:
│       • Invalid SMILES → Returns None
│       • Missing atoms → Graceful degradation
│       • Batch processing with progress tracking
│
├── 📄 toxicophores.py                     [422 lines]
│   └── Structural alert analysis
│       CLASS: ToxicophoreAnalyzer
│       20 TOXICOPHORE SMARTS PATTERNS:
│       INSECTICIDE CLASSES:
│       • Organophosphate: [P](=O)([O,S])[O,S]
│       • Carbamate: NC(=O)O
│       • Pyrethroid: CC(C)=CC(=O)
│       • Neonicotinoid: c1ncnn1 (triazole)
│       • Phenylpyrazole: n1ncc(c1)c2ccccc2
│       FUNCTIONAL GROUPS:
│       • Nitro: [N+](=O)[O-]
│       • Cyano: C#N
│       • Aromatic halogen: c[F,Cl,Br,I]
│       • Fluorinated alkyl: C(F)(F)F
│       • Sulfonyl: S(=O)(=O)
│       RING SYSTEMS:
│       • Triazole: c1ncnn1
│       • Imidazole: c1nccn1
│       • Pyridine: c1ccncc1
│       REACTIVE GROUPS:
│       • Phosphate ester: OP(=O)(O)O
│       • Aromatic amine: c-N
│       • Phenol: c[OH]
│       SPECIFIC MARKERS:
│       • Urea: NC(=O)N
│       • Oxime: C=N-O
│       • Methylenedioxyphenyl: c1cc2OCOc2cc1
│       METHODS:
│       • find_toxicophores(smiles) → Dict of matches
│       • analyze_dataset(df) → Statistical analysis
│       • calculate_enrichment(toxicophore, toxicity) → Chi-square, p-value
│       • plot_enrichment() → Bar chart
│       • save_results(path) → JSON/CSV export
│       STATISTICAL TESTS:
│       • Chi-square independence test
│       • Enrichment ratio: (toxic_with / total_with) / (toxic_without / total_without)
│       • Odds ratio calculation
│       • p-value significance (α = 0.05)
│
├── 📄 recommendations.py                  [336 lines]
│   └── KNN-based compound recommendations
│       CLASS: CompoundRecommender
│       ALGORITHM: K-Nearest Neighbors
│       WORKFLOW:
│       1. Separate toxic vs safe compounds
│       2. Fit StandardScaler on safe compounds
│       3. Build KNN index (n_neighbors=10)
│       4. For toxic compound: find K nearest safe neighbors
│       5. Rank by Euclidean distance
│       6. Convert distance → similarity: 1 / (1 + distance)
│       METHODS:
│       • fit(X_safe, X_toxic) → Train KNN
│       • recommend(compound_id, n=5) → Get alternatives
│       • batch_recommend(cid_list) → Multiple compounds
│       • calculate_similarity(c1, c2) → Pairwise similarity
│       • save_recommendations(path) → CSV export
│       USE CASES:
│       • Find safer alternatives to toxic pesticides
│       • Maintain similar molecular properties
│       • Support green chemistry initiatives
│       • Regulatory read-across
│       OUTPUT:
│       • Ranked list of safe alternatives
│       • Similarity scores (0-1)
│       • Euclidean distances
│       • Compound names and CIDs
│
├── 📄 interpretability.py                 [388 lines]
│   └── Model explainability (SHAP + LIME)
│       CLASS: ModelInterpreter
│       SHAP INTEGRATION:
│       • TreeExplainer (fast, exact for tree models)
│       • KernelExplainer (model-agnostic fallback)
│       • Global importance: Mean |SHAP|
│       • Local explanations: Waterfall plots
│       LIME INTEGRATION:
│       • TabularExplainer
│       • Per-prediction explanations
│       • Feature contribution breakdown
│       METHODS:
│       • setup_shap(explainer_type) → Initialize explainer
│       • calculate_shap_values(X) → Compute SHAP matrix
│       • plot_shap_summary() → Beeswarm plot
│       • plot_shap_importance() → Bar chart (mean |SHAP|)
│       • plot_shap_waterfall(idx) → Individual prediction
│       • explain_with_lime(idx) → LIME explanation
│       • calculate_feature_importance() → Ranking
│       • save_importance(path) → CSV export
│       VISUALIZATIONS:
│       • Summary plot (all features, all samples)
│       • Importance bar chart
│       • Waterfall plots (individual predictions)
│       • LIME HTML reports
│       TOP FEATURES IDENTIFIED:
│       1. insecticide (1.366)
│       2. herbicide (1.054)
│       3. fungicide (0.740)
│       4. year (0.641)
│       5. LogP (0.474)
│
├── 📄 temporal_analysis.py                [354 lines]
│   └── Time trend analysis
│       CLASS: TemporalAnalyzer
│       STATISTICAL TESTS:
│       • Mann-Kendall trend test (non-parametric)
│       • Sen's slope estimator
│       • Kendall's tau correlation
│       METHODS:
│       • analyze_trends(df) → Time series analysis
│       • decade_comparison() → Group by decade
│       • rolling_average(window=10) → Smoothing
│       • plot_temporal_trends() → Line charts
│       • statistical_tests() → Mann-Kendall, p-values
│       • save_results(path) → JSON export
│       INSIGHTS:
│       • Toxicity rates over 191 years (1832-2023)
│       • Decade-by-decade comparison
│       • Trend significance testing
│       • Seasonal patterns (if applicable)
│       VISUALIZATIONS:
│       • Time series plot (year vs toxicity rate)
│       • Decade boxplots
│       • Rolling average overlay
│       • Trend line with confidence intervals
│
├── 📄 chemical_space.py                   [378 lines]
│   └── Dimensionality reduction and visualization
│       CLASS: ChemicalSpaceVisualizer
│       ALGORITHMS:
│       • PCA (Principal Component Analysis)
│       • t-SNE (t-Distributed Stochastic Neighbor Embedding)
│       METHODS:
│       • fit_pca(n_components=2) → PCA transformation
│       • fit_tsne(perplexity=30) → t-SNE embedding
│       • plot_2d(color_by='toxicity') → 2D scatter
│       • plot_3d(color_by='year') → 3D interactive (Plotly)
│       • calculate_variance_explained() → PCA variance
│       • identify_clusters(method='kmeans') → Clustering
│       • save_embeddings(path) → CSV export
│       VISUALIZATIONS:
│       • PCA 2D (PC1 vs PC2, colored by toxicity)
│       • PCA 3D (PC1, PC2, PC3, interactive)
│       • t-SNE 2D (perplexity=30)
│       • Variance explained plot (scree plot)
│       INSIGHTS:
│       • Toxic vs non-toxic separation
│       • Chemical class clustering
│       • Temporal patterns in chemical space
│       • Outlier detection
│
├── 📄 source_comparison.py                [304 lines]
│   └── ECOTOX vs PPDB data source comparison
│       CLASS: SourceComparator
│       STATISTICAL TESTS:
│       • Two-sample t-test (continuous features)
│       • Chi-square test (categorical features)
│       • Mann-Whitney U test (non-parametric)
│       METHODS:
│       • compare_distributions(feature) → Statistical test
│       • plot_distribution_comparison() → Side-by-side histograms
│       • calculate_effect_size() → Cohen's d
│       • test_toxicity_agreement() → Inter-rater reliability
│       • save_comparison(path) → JSON/CSV export
│       INSIGHTS:
│       • Feature distribution differences
│       • Toxicity label agreement
│       • Sample size comparison
│       • Bias detection
│       VISUALIZATIONS:
│       • Distribution overlays (ECOTOX vs PPDB)
│       • Toxicity rate comparison
│       • Feature correlation matrices
│       • Venn diagrams (shared compounds)
│
└── 📄 __init__.py                        [Empty]
    └── Python package marker
```

---

### 📁 data/ - Datasets

```
data/
└── raw/
    └── 📄 dataset_with_descriptors.csv    [1,035 rows × 28 cols] ⭐ CORE DATA
        └── Main ApisTox dataset
            DIMENSIONS: 1,035 compounds × 28 columns
            TIME SPAN: 1832-2023 (191 years)
            SOURCES: ECOTOX (EPA), PPDB (UK)

            COLUMN SCHEMA:
            ┌────────────────────────────────────────────────────────┐
            │ IDENTIFIERS (4 columns)                                │
            ├────────────────────────────────────────────────────────┤
            │ • CID: Compound ID (integer, unique)                   │
            │ • Preferred_name: Chemical name (string)               │
            │ • SMILES: Molecular structure (string)                 │
            │ • InChI: Chemical identifier (string)                  │
            └────────────────────────────────────────────────────────┘

            ┌────────────────────────────────────────────────────────┐
            │ METADATA (4 columns)                                   │
            ├────────────────────────────────────────────────────────┤
            │ • source: ECOTOX or PPDB                               │
            │ • year: Registration/test year (1832-2023)             │
            │ • toxicity_type: Contact or Oral                       │
            │ • chemical_type: Insecticide, Herbicide, etc.          │
            └────────────────────────────────────────────────────────┘

            ┌────────────────────────────────────────────────────────┐
            │ CHEMICAL TYPE FLAGS (4 binary columns)                 │
            ├────────────────────────────────────────────────────────┤
            │ • insecticide: 1 or 0                                  │
            │ • herbicide: 1 or 0                                    │
            │ • fungicide: 1 or 0                                    │
            │ • other_agrochemical: 1 or 0                           │
            └────────────────────────────────────────────────────────┘

            ┌────────────────────────────────────────────────────────┐
            │ MOLECULAR DESCRIPTORS (15 columns, RDKit-derived)      │
            ├────────────────────────────────────────────────────────┤
            │ • MolecularWeight: 0-2000 Da                           │
            │ • LogP: -10 to 20                                      │
            │ • NumHDonors: 0-50                                     │
            │ • NumHAcceptors: 0-50                                  │
            │ • NumRotatableBonds: 0-100                             │
            │ • NumAromaticRings: 0-20                               │
            │ • TPSA: 0-500 Ų                                        │
            │ • NumHeteroatoms: 0-100                                │
            │ • NumRings: 0-20                                       │
            │ • NumSaturatedRings: 0-20                              │
            │ • NumAliphaticRings: 0-20                              │
            │ • FractionCSP3: 0-1                                    │
            │ • MolarRefractivity: 0-500                             │
            │ • BertzCT: 0-10000                                     │
            │ • HeavyAtomCount: 0-200                                │
            └────────────────────────────────────────────────────────┘

            ┌────────────────────────────────────────────────────────┐
            │ TARGET VARIABLE (1 binary column)                      │
            ├────────────────────────────────────────────────────────┤
            │ • EPA_binary: 0 (non-toxic) or 1 (toxic to bees)       │
            │   - Class 0: 739 compounds (71.4%)                     │
            │   - Class 1: 296 compounds (28.6%)                     │
            │   - Imbalance ratio: 2.5:1                             │
            └────────────────────────────────────────────────────────┘

            DATA QUALITY:
            • Missing values: <1% (imputed)
            • Duplicates: None (CID unique)
            • Outliers: Validated by domain experts
            • SMILES validity: 100% (RDKit-parseable)
```

---

### 📁 outputs/ - Generated Artifacts

```
outputs/
├── models/                                # Trained ML models
│   ├── 📄 best_model_xgboost.pkl         [13.2 MB] ⭐ PRODUCTION MODEL
│   │   └── XGBoost classifier (best performer)
│   │       • Accuracy: 83.6%
│   │       • ROC-AUC: 85.8%
│   │       • Trained on 1,035 compounds
│   │       • Hyperparameter-tuned (GridSearchCV)
│   │       • Persisted via Joblib
│   │
│   ├── 📄 best_model_random_forest.pkl   [8.7 MB]
│   │   └── Random Forest classifier (alternative)
│   │       • Accuracy: 83.2%
│   │       • 100 trees, max_depth=20
│   │
│   ├── 📄 best_model_lightgbm.pkl        [2.1 MB]
│   │   └── LightGBM classifier
│   │       • Accuracy: 82.8%
│   │       • Faster inference than XGBoost
│   │
│   ├── 📄 best_model_svm.pkl             [5.3 MB]
│   │   └── SVM classifier
│   │       • Accuracy: 81.5%
│   │       • RBF kernel
│   │
│   ├── 📄 best_model_mlp.pkl             [1.8 MB]
│   │   └── Neural network classifier
│   │       • Accuracy: 80.3%
│   │       • 3 hidden layers (100, 50, 25)
│   │
│   └── 📄 best_model_logistic.pkl        [0.3 MB]
│       └── Logistic regression (baseline)
│           • Accuracy: 78.9%
│           • Fast inference
│
├── preprocessors/                         # Data preprocessing artifacts
│   └── 📄 preprocessor.pkl               [0.5 MB] ⭐ CRITICAL
│       └── DataPreprocessor instance
│           • StandardScaler (fitted on training data)
│           • One-hot encoder settings
│           • Feature names and order
│           • Class label mapping
│           • Required for all predictions
│
├── metrics/                               # Model performance metrics
│   ├── 📄 training_results.json          [8 KB]
│   │   └── Model comparison results
│   │       • All 6 models' metrics
│   │       • Accuracy, Precision, Recall, F1, ROC-AUC
│   │       • Training time, parameters
│   │       • Cross-validation scores
│   │
│   ├── 📄 feature_importance_shap.csv    [2 KB]
│   │   └── SHAP feature importance
│   │       • 24 features ranked by mean |SHAP|
│   │       • Standard deviations
│   │       • Top features:
│   │         1. insecticide (1.366)
│   │         2. herbicide (1.054)
│   │         3. fungicide (0.740)
│   │         4. year (0.641)
│   │         5. LogP (0.474)
│   │
│   ├── 📄 confusion_matrix.csv           [1 KB]
│   │   └── Test set confusion matrix
│   │       • True Positives, False Positives
│   │       • True Negatives, False Negatives
│   │
│   └── 📄 classification_report.txt      [1 KB]
│       └── Sklearn classification report
│           • Per-class precision, recall, F1
│           • Support (sample counts)
│           • Weighted and macro averages
│
├── figures/                               # Visualizations (20+ PNGs)
│   ├── 📄 shap_summary.png               [1.2 MB]
│   │   └── SHAP beeswarm plot
│   │       • All features, all samples
│   │       • Color: feature value (blue=low, red=high)
│   │       • X-axis: SHAP value (impact on prediction)
│   │
│   ├── 📄 shap_importance.png            [0.8 MB]
│   │   └── SHAP feature importance bar chart
│   │       • Mean |SHAP| values
│   │       • Sorted descending
│   │
│   ├── 📄 toxicophore_enrichment.png     [0.9 MB]
│   │   └── Toxicophore enrichment bar chart
│   │       • Enrichment ratio for each toxicophore
│   │       • Error bars (95% CI)
│   │       • Statistical significance markers
│   │
│   ├── 📄 toxicophore_prevalence.png     [0.7 MB]
│   │   └── Prevalence vs toxicity scatter
│   │       • X: prevalence (% compounds with toxicophore)
│   │       • Y: toxicity rate (% toxic among those)
│   │       • Bubble size: sample size
│   │
│   ├── 📄 chemical_space_pca.png         [1.5 MB]
│   │   └── PCA 2D scatter plot
│   │       • PC1 vs PC2
│   │       • Color: toxic (red) vs non-toxic (green)
│   │       • Variance explained labels
│   │
│   ├── 📄 chemical_space_tsne.png        [1.3 MB]
│   │   └── t-SNE 2D embedding
│   │       • Perplexity=30
│   │       • Color: toxicity
│   │       • Clustering visible
│   │
│   ├── 📄 temporal_trends.png            [1.1 MB]
│   │   └── Toxicity rate over time
│   │       • X: year (1832-2023)
│   │       • Y: toxicity rate (%)
│   │       • Trend line (Mann-Kendall)
│   │       • Rolling average overlay
│   │
│   ├── 📄 decade_comparison.png          [0.9 MB]
│   │   └── Decade boxplots
│   │       • Toxicity rate by decade
│   │       • Statistical significance markers
│   │
│   ├── 📄 correlation_matrix.png         [2.0 MB]
│   │   └── Feature correlation heatmap
│   │       • All 24 features
│   │       • Color: correlation (-1 to 1)
│   │       • Hierarchical clustering
│   │
│   ├── 📄 confusion_matrix_plot.png      [0.6 MB]
│   │   └── Confusion matrix heatmap
│   │       • Test set predictions
│   │       • Annotated with counts
│   │
│   ├── 📄 roc_curve.png                  [0.7 MB]
│   │   └── ROC curve plot
│   │       • True Positive Rate vs False Positive Rate
│   │       • AUC = 0.858
│   │       • Diagonal reference line
│   │
│   ├── 📄 precision_recall_curve.png     [0.6 MB]
│   │   └── Precision-Recall curve
│   │       • Precision vs Recall
│   │       • AP (Average Precision) score
│   │
│   └── (8 more waterfall plots for individual predictions)
│       ├── waterfall_0.png ... waterfall_7.png
│       └── SHAP waterfall plots (local explanations)
│
└── analysis/                              # Analysis results
    ├── 📄 toxicophore_results.json       [12 KB]
    │   └── Toxicophore analysis results
    │       • 20 toxicophores analyzed
    │       • Prevalence, enrichment, p-values
    │       • Chi-square test results
    │       • Odds ratios
    │       TOP TOXICOPHORES (enrichment):
    │       1. Organophosphate (3.2× enrichment, p<0.001)
    │       2. Carbamate (2.8× enrichment, p<0.001)
    │       3. Neonicotinoid (2.5× enrichment, p<0.01)
    │
    ├── 📄 alternatives.csv               [45 KB]
    │   └── KNN safer alternatives
    │       • For each toxic compound:
    │         - Top 5 safe alternatives
    │         - Similarity scores
    │         - Euclidean distances
    │         - Compound names, CIDs
    │       • 296 toxic compounds × 5 alternatives = 1,480 rows
    │
    ├── 📄 temporal_trends.json           [8 KB]
    │   └── Temporal analysis results
    │       • Mann-Kendall test results
    │       • Sen's slope estimate
    │       • p-value for trend significance
    │       • Decade-by-decade statistics
    │       • Rolling averages (10-year window)
    │
    └── 📄 chemical_space_results.json    [15 KB]
        └── Chemical space analysis
            • PCA variance explained (PC1: 23%, PC2: 18%, ...)
            • t-SNE hyperparameters
            • Cluster assignments (K-means)
            • Silhouette scores
            • Outlier detection results
```

---

### 📁 tests/ - Unit and Integration Tests

```
tests/
├── 📄 test_api.py                        [250 lines]
│   └── API endpoint tests (pytest)
│       TEST CASES (15):
│       • test_root_endpoint() → GET /
│       • test_health_check() → GET /health
│       • test_predict_valid_input() → POST /predict (valid)
│       • test_predict_invalid_input() → POST /predict (invalid)
│       • test_predict_smiles_valid() → POST /predict/smiles (valid)
│       • test_predict_smiles_invalid() → POST /predict/smiles (bad SMILES)
│       • test_model_info() → GET /model/info
│       • test_feature_importance() → GET /feature/importance
│       • test_history() → GET /history
│       • test_toxicophores() → GET /analysis/toxicophores
│       • test_toxicophores_molecule() → POST /analysis/toxicophores/molecule
│       • test_recommend_alternatives() → GET /recommend/alternatives/{cid}
│       • test_cors_headers() → CORS middleware
│       • test_error_handling_404() → 404 Not Found
│       • test_error_handling_500() → 500 Internal Server Error
│       FIXTURES:
│       • client: FastAPI TestClient
│       • valid_input: Example PredictionInput
│       • invalid_input: Malformed input (for validation testing)
│       MOCKING:
│       • Model predictions (to avoid loading actual model)
│       • Preprocessor loading
│       ASSERTIONS:
│       • Status codes (200, 400, 404, 422, 500)
│       • Response schema validation
│       • JSON structure
│
├── 📄 test_models.py                     [200 lines]
│   └── Model training tests
│       TEST CASES (10):
│       • test_model_instantiation() → get_model() for all 6
│       • test_hyperparameter_tuning() → GridSearchCV
│       • test_cross_validation() → 5-fold CV
│       • test_metrics_calculation() → Accuracy, F1, ROC-AUC
│       • test_model_persistence() → save_model() / load_model()
│       • test_prediction_consistency() → Same input → same output
│       • test_train_all_models() → Train and compare 6 models
│       • test_overfitting_detection() → Train vs validation gap
│       • test_class_imbalance_handling() → SMOTE effect
│       • test_model_comparison() → Best model selection
│       FIXTURES:
│       • sample_data: Small dataset (100 compounds)
│       • trained_model: Pre-trained XGBoost
│       ASSERTIONS:
│       • Metrics > baseline (accuracy > 60%)
│       • Model files created
│       • Loaded model == saved model
│
├── 📄 test_preprocessing.py              [180 lines]
│   └── Preprocessing pipeline tests
│       TEST CASES (12):
│       • test_data_loading() → load_data()
│       • test_feature_separation() → prepare_features()
│       • test_categorical_encoding() → encode_categorical_features()
│       • test_train_test_split() → split_data()
│       • test_stratification() → Class distribution preserved
│       • test_feature_scaling() → StandardScaler
│       • test_smote_balancing() → handle_imbalance()
│       • test_scaffold_split() → scaffold_split()
│       • test_preprocessor_persistence() → save/load
│       • test_missing_value_handling() → Imputation
│       • test_outlier_detection() → IQR method
│       • test_pipeline_integration() → End-to-end
│       FIXTURES:
│       • sample_csv: Temporary CSV file
│       • preprocessor: DataPreprocessor instance
│       ASSERTIONS:
│       • Split sizes correct (70/10/20)
│       • Class ratios preserved
│       • Scaled features: mean≈0, std≈1
│       • SMOTE: balanced classes
│
└── 📄 __init__.py                        [Empty]
    └── Python package marker
```

---

### 📁 docs/ - Documentation

```
docs/
├── 📄 API_DOCS.md                        [500 lines]
│   └── Comprehensive API documentation
│       • All 10 endpoints documented
│       • Request/response schemas
│       • Example curl commands
│       • Authentication (future)
│       • Rate limiting (future)
│       • Error codes and handling
│       • Pagination (history endpoint)
│
├── 📄 MODEL_CARD.md                      [300 lines]
│   └── Model documentation (ML best practice)
│       SECTIONS:
│       • Model Details (algorithm, version, authors)
│       • Intended Use (toxicity prediction for bees)
│       • Training Data (ApisTox dataset description)
│       • Performance Metrics (accuracy, ROC-AUC, etc.)
│       • Limitations (data imbalance, domain applicability)
│       • Ethical Considerations (regulatory use, safety)
│       • Citation (how to cite)
│
├── 📄 project_proposal.md                [250 lines]
│   └── Academic project proposal
│       • Background and motivation
│       • Research questions
│       • Methodology
│       • Expected outcomes
│       • Timeline
│       • References
│
├── 📄 ARCHITECTURE_OVERVIEW.md           [AUTO-GENERATED] ⭐ NEW
│   └── System architecture documentation
│       • High-level architecture diagrams
│       • Component interaction maps
│       • Data flow diagrams
│       • Technology stack
│       • Deployment architecture
│       • Security measures
│       • Performance characteristics
│
├── 📄 DIRECTORY_STRUCTURE.md             [THIS FILE] ⭐ NEW
│   └── Annotated directory tree
│
└── 📄 FILES_INVENTORY.md                 [PENDING] ⭐ NEW
    └── Critical files inventory (to be generated)
```

---

## File Size Summary

### By Category

| Category | File Count | Total Size | Percentage |
|----------|------------|------------|------------|
| **Models** (.pkl) | 6 | ~31 MB | 68% |
| **Figures** (.png) | 20+ | ~12 MB | 26% |
| **Data** (.csv) | 1 | ~0.8 MB | 2% |
| **Code** (.py, .tsx, .ts) | 25 | ~0.5 MB | 1% |
| **Config** (.json, .yml) | 15 | ~0.2 MB | <1% |
| **Docs** (.md) | 8 | ~0.3 MB | 1% |
| **Tests** (.py) | 3 | ~0.1 MB | <1% |
| **Other** | 10+ | ~0.5 MB | 1% |

**Total Project Size**: ~45 MB (excluding node_modules, venv, caches)

---

## File Naming Conventions

### Python Files
- **Snake case**: `molecular_features.py`, `temporal_analysis.py`
- **Descriptive**: Names indicate purpose (preprocessing, models, etc.)
- **No abbreviations**: Full words preferred

### TypeScript/React Files
- **PascalCase** (components): `PredictionForm.tsx`, `ResultDisplay.tsx`
- **camelCase** (utilities): `api.ts`
- **Extension**: `.tsx` for JSX, `.ts` for pure TypeScript

### Data Files
- **Snake case**: `dataset_with_descriptors.csv`
- **Descriptive**: Indicates content and format
- **Version suffixes** (if needed): `data_v2.csv`

### Configuration Files
- **Lowercase with hyphens**: `docker-compose.yml`
- **Dots for scopes**: `tsconfig.json`, `vite.config.ts`
- **Standard names**: `requirements.txt`, `package.json` (ecosystem conventions)

### Generated Artifacts
- **Snake case**: `best_model_xgboost.pkl`, `shap_summary.png`
- **Descriptive prefixes**: `best_` (models), `shap_` (SHAP plots)
- **Sequential numbering**: `waterfall_0.png`, `waterfall_1.png`, ...

---

## Critical Path Files

**For Running API**:
1. `app/backend/main.py` - API entry point
2. `outputs/models/best_model_xgboost.pkl` - Production model
3. `outputs/preprocessors/preprocessor.pkl` - Preprocessor
4. `requirements-production.txt` - Dependencies

**For Running Frontend**:
1. `app/frontend/src/main.tsx` - Entry point
2. `app/frontend/src/App.tsx` - Root component
3. `app/frontend/package.json` - Dependencies
4. `app/frontend/vite.config.ts` - Build config

**For Training Models**:
1. `data/raw/dataset_with_descriptors.csv` - Dataset
2. `src/preprocessing.py` - Preprocessing pipeline
3. `src/models.py` - Training script
4. `requirements.txt` - All dependencies

**For Deployment**:
1. `docker-compose.yml` - Container orchestration
2. `Dockerfile.backend` - Backend container
3. `Dockerfile.frontend` - Frontend container
4. `vercel.json` - Serverless config

---

**Document Maintained By**: ApisTox Development Team
**Last Review**: November 19, 2025
**Next Review**: December 2025
