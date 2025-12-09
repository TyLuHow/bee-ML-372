# Use official Python runtime
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc g++ libgomp1 libxrender1 \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python dependencies
COPY backend/requirements.txt ./backend/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r backend/requirements.txt

# Copy application
COPY backend ./backend

# Set working directory
WORKDIR /app/backend

# Expose port
EXPOSE 8080

# Start uvicorn directly (no script, no cd, no exec)
CMD python -m uvicorn api.main:app --host 0.0.0.0 --port ${PORT:-8080}
