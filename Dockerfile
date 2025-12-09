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
