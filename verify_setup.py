#!/usr/bin/env python3
"""
Quick Setup Verification Script
================================
Checks if all components are properly installed and configured.

Usage: python verify_setup.py
"""

import sys
import os
from pathlib import Path

def check_python_version():
    """Check Python version."""
    version = sys.version_info
    print(f"✓ Python {version.major}.{version.minor}.{version.micro}")
    if version.major < 3 or (version.major == 3 and version.minor < 9):
        print("  ⚠️  Warning: Python 3.9+ recommended")
        return False
    return True

def check_dependencies():
    """Check if required packages are installed."""
    required = [
        'pandas', 'numpy', 'scikit-learn', 'xgboost', 
        'fastapi', 'uvicorn', 'rdkit', 'matplotlib', 
        'seaborn', 'joblib', 'scipy', 'plotly'
    ]
    
    missing = []
    for package in required:
        try:
            __import__(package.replace('-', '_'))
            print(f"  ✓ {package}")
        except ImportError:
            print(f"  ✗ {package} - MISSING")
            missing.append(package)
    
    return len(missing) == 0, missing

def check_file_structure():
    """Check if required files and directories exist."""
    required_files = [
        'outputs/dataset_final.csv',
        'src/preprocessing.py',
        'src/models.py',
        'src/molecular_features.py',
        'app/backend/main.py',
        'requirements.txt'
    ]
    
    required_dirs = [
        'outputs',
        'outputs/models',
        'outputs/figures',
        'outputs/analysis',
        'src',
        'app/backend',
        'tests'
    ]
    
    all_ok = True
    
    print("\n📁 Checking file structure...")
    for fpath in required_files:
        if Path(fpath).exists():
            print(f"  ✓ {fpath}")
        else:
            print(f"  ✗ {fpath} - MISSING")
            all_ok = False
    
    for dpath in required_dirs:
        if Path(dpath).exists():
            print(f"  ✓ {dpath}/")
        else:
            print(f"  ⚠️  {dpath}/ - Will be created")
            os.makedirs(dpath, exist_ok=True)
    
    return all_ok

def check_models():
    """Check if models are trained."""
    print("\n🤖 Checking trained models...")
    
    model_files = [
        'outputs/models/best_model.pkl',
        'outputs/models/preprocessor.pkl'
    ]
    
    trained = True
    for mpath in model_files:
        if Path(mpath).exists():
            size = Path(mpath).stat().st_size / 1024  # KB
            print(f"  ✓ {mpath} ({size:.1f} KB)")
        else:
            print(f"  ✗ {mpath} - NOT TRAINED")
            trained = False
    
    if not trained:
        print("\n  ℹ️  Run training: python src/models.py")
    
    return trained

def check_analyses():
    """Check if analyses have been run."""
    print("\n📊 Checking analysis outputs...")
    
    analysis_files = [
        'outputs/figures/temporal_trends.png',
        'outputs/figures/chemical_space_pca.png',
        'outputs/analysis/temporal_analysis.json'
    ]
    
    completed = 0
    for apath in analysis_files:
        if Path(apath).exists():
            print(f"  ✓ {apath}")
            completed += 1
        else:
            print(f"  ✗ {apath} - NOT GENERATED")
    
    if completed < len(analysis_files):
        print("\n  ℹ️  Run analyses: python run_comprehensive_analysis.py")
    
    return completed == len(analysis_files)

def check_api():
    """Check if API can be imported."""
    print("\n🌐 Checking API...")
    try:
        from app.backend.main import app
        print("  ✓ API module imports successfully")
        return True
    except Exception as e:
        print(f"  ✗ API import failed: {e}")
        return False

def main():
    print("="*60)
    print("🐝 APISTOX SETUP VERIFICATION")
    print("="*60)
    
    print("\n🐍 Checking Python environment...")
    py_ok = check_python_version()
    
    print("\n📦 Checking dependencies...")
    deps_ok, missing = check_dependencies()
    
    file_ok = check_file_structure()
    model_ok = check_models()
    analysis_ok = check_analyses()
    api_ok = check_api()
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    checks = {
        "Python version": py_ok,
        "Dependencies": deps_ok,
        "File structure": file_ok,
        "Trained models": model_ok,
        "Analysis outputs": analysis_ok,
        "API module": api_ok
    }
    
    for check, status in checks.items():
        symbol = "✅" if status else "❌"
        print(f"{symbol} {check}")
    
    all_ok = all(checks.values())
    
    if not all_ok:
        print("\n" + "="*60)
        print("⚠️  SETUP INCOMPLETE - Follow these steps:")
        print("="*60)
        
        if not deps_ok:
            print("\n1. Install missing dependencies:")
            print(f"   pip install {' '.join(missing)}")
        
        if not model_ok:
            print("\n2. Train models:")
            print("   python src/models.py")
        
        if not analysis_ok:
            print("\n3. Run analyses:")
            print("   python run_comprehensive_analysis.py")
        
        print("\n📖 See SETUP_GUIDE.md for detailed instructions")
    else:
        print("\n" + "="*60)
        print("✅ SETUP COMPLETE!")
        print("="*60)
        print("\n🚀 Next steps:")
        print("   1. Start API: python -m uvicorn app.backend.main:app --reload --port 8000")
        print("   2. Start frontend: cd app/frontend && npm run dev")
        print("   3. Open browser: http://localhost:5173")
        print("\n📚 Documentation: http://localhost:8000/docs")
    
    return 0 if all_ok else 1

if __name__ == "__main__":
    sys.exit(main())




