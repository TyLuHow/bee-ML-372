@echo off
REM Quick Start Script for ApisTox
REM ================================
REM This script automates the initial setup process
REM Run this ONCE after cloning the repository

echo =====================================
echo APISTOX QUICK START
echo =====================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found. Please install Python 3.9+ first.
    pause
    exit /b 1
)

echo Step 1: Installing Python dependencies...
echo.
pip install -r requirements.txt
if errorlevel 1 (
    echo ERROR: Failed to install dependencies
    pause
    exit /b 1
)
echo ✓ Dependencies installed
echo.

echo Step 2: Verifying setup...
echo.
python verify_setup.py
echo.

echo Step 3: Training models (this may take 5-10 minutes)...
echo.
python src\models.py
if errorlevel 1 (
    echo ERROR: Model training failed
    pause
    exit /b 1
)
echo ✓ Models trained
echo.

echo Step 4: Running comprehensive analysis...
echo.
python run_comprehensive_analysis.py
if errorlevel 1 (
    echo WARNING: Some analyses may have failed, but continuing...
)
echo ✓ Analyses complete
echo.

echo =====================================
echo ✅ SETUP COMPLETE!
echo =====================================
echo.
echo Next steps:
echo   1. Start the API:
echo      python -m uvicorn app.backend.main:app --reload --port 8000
echo.
echo   2. In a NEW terminal, start the frontend:
echo      cd app\frontend
echo      npm install
echo      npm run dev
echo.
echo   3. Open your browser to: http://localhost:5173
echo.
echo See SETUP_GUIDE.md for more information
echo.
pause




