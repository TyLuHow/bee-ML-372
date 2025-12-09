<objective>
Conduct a comprehensive audit of the ApisTox deployment infrastructure and machine learning model integration. This analysis will identify deployment issues, understand the current state of services across Vercel, Railway, and GitHub, and provide a detailed understanding of how the ML model is integrated and used throughout the application.
</objective>

<context>
This is the ApisTox (bee-ml-372) project - a bee toxicity prediction platform with:
- Frontend: React/TypeScript deployed on Vercel
- Backend: Python ML service deployed on Railway
- Repository: GitHub at TyLuHow/bee-ml-372

The analysis will inform decisions about fixing broken deployments and understanding the ML architecture for future improvements.
</context>

<phase_1_deployment_audit>
<title>Multi-Platform Deployment Status Check</title>

Thoroughly analyze the health and configuration of all deployment platforms:

1. **Vercel Analysis**
   - Check deployment status: `vercel list --project bee-ml-372`
   - Get latest deployment details: `vercel inspect [deployment-url]`
   - Review build logs for errors or warnings
   - Verify environment variables are properly configured
   - Check domain configurations and SSL status
   - Identify any failed builds or deployments

2. **Railway Analysis**
   - Check service status: `railway status`
   - List all services: `railway list`
   - Get deployment logs: `railway logs` (look for errors, crashes, startup issues)
   - Review railway.json and railway.toml configurations
   - Verify environment variables and secrets
   - Check health endpoints if available
   - Identify resource usage and any scaling issues

3. **GitHub Repository Analysis**
   - Check recent workflow runs: `gh run list --repo TyLuHow/bee-ml-372`
   - Review failed workflows: `gh run view [run-id]` for any failures
   - Check branch status and protection rules
   - Review recent commits for deployment-related changes
   - Identify any CI/CD configuration issues

4. **Cross-Platform Integration**
   - Verify frontend can connect to backend API
   - Check CORS configurations
   - Test API endpoints from frontend deployment
   - Identify any environment-specific issues (prod vs preview)
</phase_1_deployment_audit>

<phase_2_model_integration_analysis>
<title>ML Model Usage and Architecture Deep Dive</title>

Thoroughly explore and document how the machine learning model is used in the application:

1. **Model Files and Artifacts**
   - Locate model files (examine @backend for .pkl, .h5, .pt, .onnx, or similar)
   - Identify model type and framework (scikit-learn, TensorFlow, PyTorch, etc.)
   - Check model version and last updated date
   - Document model size and storage location

2. **Model Loading and Initialization**
   - Find where model is loaded in the backend code (examine @backend/\*\*/\*.py)
   - Identify initialization logic and startup procedures
   - Check for model caching or lazy loading strategies
   - Review error handling for model loading failures

3. **Inference Pipeline**
   - Trace the complete prediction flow from API request to response
   - Identify preprocessing steps (data cleaning, feature engineering)
   - Document the actual inference/prediction call
   - Review postprocessing steps (result formatting, confidence scores)
   - Map out all files involved in the prediction pipeline

4. **API Endpoints**
   - List all model-related API endpoints (examine @backend for route definitions)
   - Document request/response schemas
   - Check rate limiting or caching strategies
   - Identify authentication/authorization for model access

5. **Frontend Integration**
   - Find where frontend calls the model API (examine @App.tsx and @services/\*)
   - Trace user input flow to prediction request
   - Review how predictions are displayed to users
   - Check error handling and loading states

6. **Model Metadata and Configuration**
   - Look for model configuration files (examine @backend for config.py, settings.py, etc.)
   - Document hyperparameters or model settings
   - Check feature names and expected input schema
   - Identify any A/B testing or model versioning logic

7. **Performance and Monitoring**
   - Check for prediction latency logging
   - Look for model performance metrics collection
   - Identify any monitoring or alerting for model health
</phase_2_model_integration_analysis>

<analysis_approach>
For maximum efficiency, whenever you need to perform multiple independent operations (checking different platforms, reading multiple files), invoke all relevant tools simultaneously rather than sequentially.

After receiving tool results, carefully reflect on their quality and determine optimal next steps before proceeding. If you encounter authentication issues with CLI tools, document them and suggest manual verification steps.
</analysis_approach>

<output>
Create a comprehensive analysis report at: `./analysis/deployment-and-model-audit.md`

The report should include:

1. **Executive Summary**
   - Overall health status (working/broken/degraded)
   - Critical issues requiring immediate attention
   - Quick wins for improvement

2. **Deployment Status by Platform**
   - Vercel: status, issues, configuration notes
   - Railway: status, issues, configuration notes
   - GitHub: workflow status, recent activity
   - Cross-platform connectivity and integration issues

3. **ML Model Architecture**
   - Model type, framework, and specifications
   - Complete data flow diagram (user input → preprocessing → inference → response)
   - File-by-file breakdown of model usage with file paths and line numbers
   - API endpoint documentation
   - Frontend integration points

4. **Issues and Recommendations**
   - Deployment issues with severity (critical/high/medium/low)
   - Model integration concerns or optimization opportunities
   - Security considerations
   - Performance improvement suggestions

5. **Action Items**
   - Prioritized list of fixes needed
   - Steps to resolve each issue
   - Long-term improvements to consider
</output>

<verification>
Before completing, verify your analysis:
- All three platforms (Vercel, Railway, GitHub) have been checked
- Model files have been located and examined
- Complete prediction flow has been traced from frontend to backend
- All findings are documented with specific file paths and line numbers
- Recommendations are actionable and prioritized
</verification>

<success_criteria>
- Comprehensive deployment status for all platforms with specific error messages if issues exist
- Complete understanding of ML model architecture with code references
- Clear action items for fixing any broken deployments
- Documented model usage pattern that could guide future development
</success_criteria>
