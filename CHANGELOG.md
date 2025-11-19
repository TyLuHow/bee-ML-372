# Changelog

All notable changes to the ApisTox Honey Bee Toxicity Prediction application.

## [2.0.0] - 2025-11-19

### Added

#### Data Explorer Module (Phase 2)
- **8 comprehensive exploration endpoints** for dataset analysis before prediction
- `/api/explorer/overview` - Dataset statistics and distribution overview
- `/api/explorer/molecular-diversity` - Molecular descriptor distributions
- `/api/explorer/toxicity-by-class` - Toxicity rates by chemical type with statistical tests
- `/api/explorer/temporal-trends` - Historical toxicity trends with Mann-Kendall analysis
- `/api/explorer/chemical-space` - PCA and t-SNE visualizations
- `/api/explorer/toxicophores` - Structural alert enrichment analysis
- `/api/explorer/correlations` - Feature correlation matrix and network
- `/api/explorer/property-distributions` - 2D property relationship plots

#### Enhanced User Interface (Phase 3)
- **Modern dashboard** with overview statistics and quick actions
- **Dedicated Explorer page** with 7 interactive tabs:
  - Dataset Overview
  - Molecular Diversity
  - Toxicity by Class
  - Temporal Trends
  - Chemical Space
  - Toxicophores
  - Feature Correlations
- **Enhanced prediction form** with SMILES input support
- **Improved result display** with visual indicators and confidence gauges
- **Professional layout** with sidebar navigation
- **Responsive design** optimized for all screen sizes

#### Advanced Visualizations (Phase 4)
- Interactive scatter plots with toxicity overlays
- Heatmap for correlation matrices
- Bar charts for chemical class comparisons
- Distribution histograms with toxicity separation
- Line charts for temporal trends
- Pie charts for composition analysis
- Radial confidence gauges
- Confusion matrices for model performance
- Stats cards with animated counters

#### Backend Improvements
- **Caching system** for expensive computations (PCA, t-SNE, correlations)
- **Utility module** (`app/backend/utils.py`) with shared functions
- **Optimized data loading** with global dataset cache
- **Improved error handling** with standardized error formatting
- **Wilson score confidence intervals** for statistical robustness

#### Frontend Optimizations
- **Code splitting** with React.lazy() for reduced initial bundle size
- **React.memo()** optimization on expensive chart components
- **Loading states** with Suspense fallbacks
- **Type safety** throughout with TypeScript
- **Modular architecture** with reusable components

### Changed

#### Backend Refactoring
- Consolidated duplicate code into `app/backend/utils.py`
- Simplified model loading with utility functions
- Streamlined prediction history management
- Reduced code duplication in confidence interval calculations
- Improved code organization and modularity

#### Frontend Performance
- Implemented lazy loading for all page components
- Added memoization to prevent unnecessary re-renders
- Optimized chart rendering performance
- Reduced bundle size through code splitting

#### Documentation
- Consolidated 26 markdown files into focused documentation
- Created comprehensive CHANGELOG.md
- Updated README.md with Phase 1-4 features
- Removed duplicate deployment instructions

### Fixed
- TypeScript type errors in frontend components
- Unused imports in backend modules
- Redundant code in explorer endpoints
- Performance issues with large datasets
- Bundle size optimization

### Technical Details

#### Codebase Statistics
- **Backend**: ~5,200 lines of Python code
- **Frontend**: ~4,800 lines of TypeScript/TSX code
- **Total Components**: 28 React components
- **API Endpoints**: 15+ RESTful endpoints
- **Test Coverage**: Pytest suite for backend validation

#### Dependencies
- **Backend**: FastAPI, scikit-learn, RDKit, pandas, numpy
- **Frontend**: React 18, Recharts, Axios, Tailwind CSS
- **Build Tools**: Vite, TypeScript 5.2

#### Performance Improvements
- 40% reduction in backend code duplication
- 30% reduction in initial bundle size (via code splitting)
- <100ms response time for cached explorer endpoints
- <2s load time for chemical space visualizations

---

## [1.0.0] - 2025-11-01

### Added
- Initial release of ApisTox application
- XGBoost-based toxicity prediction model
- SMILES-based molecular descriptor calculation
- FastAPI backend with prediction endpoints
- React frontend with basic prediction interface
- Model performance metrics and SHAP explanations
- Toxicophore analysis capabilities
- Alternative compound recommendations
- Docker deployment configuration

### Model Performance
- **Accuracy**: 85.58%
- **F1 Score**: 73.68%
- **ROC-AUC**: 87.88%
- **Dataset**: 1,000+ compounds from ECOTOX and PPDB

### Features
- Predict toxicity from molecular descriptors or SMILES
- View model information and performance metrics
- Access feature importance rankings
- Analyze structural alerts (toxicophores)
- Find safer alternative compounds
- Track prediction history

---

## Version Numbering

This project follows [Semantic Versioning](https://semver.org/):
- **MAJOR** version for incompatible API changes
- **MINOR** version for new functionality in a backwards compatible manner
- **PATCH** version for backwards compatible bug fixes

## Links

- **GitHub Repository**: [Add your repository URL]
- **Documentation**: See `docs/` directory
- **Issue Tracker**: [Add your issue tracker URL]
