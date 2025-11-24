#!/usr/bin/env python3
"""
System Integration Test Script
===============================

Tests the complete system workflow:
1. Model and preprocessor loading
2. Data preprocessing
3. Model prediction
4. API functionality (if running)

Author: IME 372 Project Team
Date: November 2025
"""

import os
import sys
import joblib
import pandas as pd
import numpy as np
from pathlib import Path

def test_model_files():
    """Test that model and preprocessor files exist and can be loaded."""
    print("=" * 80)
    print("TEST 1: Model and Preprocessor Files")
    print("=" * 80)
    
    model_path = "outputs/models/best_model_xgboost.pkl"
    prep_path = "outputs/preprocessors/preprocessor.pkl"
    
    # Check existence
    model_exists = os.path.exists(model_path)
    prep_exists = os.path.exists(prep_path)
    
    print(f"Model file exists: {model_exists}")
    print(f"Preprocessor file exists: {prep_exists}")
    
    if not model_exists:
        print("[FAIL] Model file not found!")
        return False
    
    if not prep_exists:
        print("[FAIL] Preprocessor file not found!")
        return False
    
    # Try loading
    try:
        model = joblib.load(model_path)
        print(f"[OK] Model loaded successfully: {type(model).__name__}")
    except Exception as e:
        print(f"[FAIL] Error loading model: {e}")
        return False
    
    try:
        # Fix pickle module path issue (preprocessing -> src.preprocessing)
        import sys
        import src.preprocessing as preprocessing_module
        sys.modules['preprocessing'] = preprocessing_module
        
        prep = joblib.load(prep_path)
        print(f"[OK] Preprocessor loaded successfully: {type(prep).__name__}")
    except Exception as e:
        print(f"[FAIL] Error loading preprocessor: {e}")
        return False
    
    print("[PASS] Test 1 PASSED\n")
    return True


def test_data_loading():
    """Test that dataset can be loaded."""
    print("=" * 80)
    print("TEST 2: Dataset Loading")
    print("=" * 80)
    
    dataset_path = "outputs/dataset_final.csv"
    
    if not os.path.exists(dataset_path):
        print(f"[FAIL] Dataset not found: {dataset_path}")
        return False
    
    try:
        df = pd.read_csv(dataset_path)
        print(f"Dataset loaded: {df.shape[0]} rows × {df.shape[1]} columns")
        print(f"Columns: {list(df.columns)}")
        print(f"Target distribution:\n{df['label'].value_counts()}")
        print("[PASS] Test 2 PASSED\n")
        return True
    except Exception as e:
        print(f"[FAIL] Error loading dataset: {e}")
        return False


def test_prediction_pipeline():
    """Test complete prediction pipeline."""
    print("=" * 80)
    print("TEST 3: Prediction Pipeline")
    print("=" * 80)
    
    try:
        # Fix pickle module path issue (preprocessing -> src.preprocessing)
        import sys
        import src.preprocessing as preprocessing_module
        sys.modules['preprocessing'] = preprocessing_module

        # Load model and preprocessor
        model = joblib.load("outputs/models/best_model_xgboost.pkl")
        preprocessor = joblib.load("outputs/preprocessors/preprocessor.pkl")
        
        # Create sample input
        sample_input = pd.DataFrame({
            'source': ['PPDB'],
            'year': [2020],
            'toxicity_type': ['Contact'],
            'herbicide': [0],
            'fungicide': [0],
            'insecticide': [1],
            'other_agrochemical': [0],
            'MolecularWeight': [350.5],
            'LogP': [3.2],
            'NumHDonors': [2],
            'NumHAcceptors': [4],
            'NumRotatableBonds': [5],
            'AromaticRings': [2],
            'TPSA': [65.3],
            'NumHeteroatoms': [5],
            'NumAromaticAtoms': [12],
            'NumSaturatedRings': [0],
            'NumAliphaticRings': [0],
            'RingCount': [2],
            'FractionCsp3': [0.25],
            'NumAromaticCarbocycles': [1],
            'NumSaturatedCarbocycles': [0]
        })
        
        print("Sample input created:")
        print(sample_input.T)
        
        # Preprocess
        X_processed = preprocessor.transform(sample_input)
        print(f"\nProcessed input shape: {X_processed.shape}")
        
        # Predict
        prediction = model.predict(X_processed)
        probability = model.predict_proba(X_processed)
        
        print(f"\nPrediction: {prediction[0]} ({'Toxic' if prediction[0] == 1 else 'Non-Toxic'})")
        print(f"Probability [Non-Toxic, Toxic]: {probability[0]}")
        print(f"Confidence: {max(probability[0]):.2%}")
        
        print("[PASS] Test 3 PASSED\n")
        return True
        
    except Exception as e:
        print(f"[FAIL] Error in prediction pipeline: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_api_imports():
    """Test that API dependencies are available."""
    print("=" * 80)
    print("TEST 4: API Dependencies")
    print("=" * 80)
    
    try:
        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from pydantic import BaseModel
        print("[OK] FastAPI imported")
        print("[OK] TestClient imported")
        print("[OK] Pydantic imported")
        print("[PASS] Test 4 PASSED\n")
        return True
    except ImportError as e:
        print(f"[FAIL] Missing API dependency: {e}")
        print("Install with: pip install fastapi pydantic python-multipart")
        return False


def test_visualization_libraries():
    """Test that visualization libraries are available."""
    print("=" * 80)
    print("TEST 5: Visualization Libraries")
    print("=" * 80)
    
    try:
        import matplotlib
        print(f"[OK] matplotlib {matplotlib.__version__}")
        import seaborn
        print(f"[OK] seaborn {seaborn.__version__}")
        print("[PASS] Test 5 PASSED\n")
        return True
    except ImportError as e:
        print(f"[WARN] Visualization library missing: {e}")
        print("Install with: pip install matplotlib seaborn")
        return False


def test_interpretability_libraries():
    """Test that interpretability libraries are available."""
    print("=" * 80)
    print("TEST 6: Interpretability Libraries")
    print("=" * 80)
    
    has_shap = False
    has_lime = False
    
    try:
        import shap
        print(f"[OK] SHAP {shap.__version__}")
        has_shap = True
    except ImportError:
        print("[WARN] SHAP not installed (pip install shap)")
    
    try:
        import lime
        print(f"[OK] LIME")
        has_lime = True
    except ImportError:
        print("[WARN] LIME not installed (pip install lime)")
    
    if has_shap or has_lime:
        print("[PASS] Test 6 PASSED (partial)\n")
        return True
    else:
        print("[FAIL] Test 6 FAILED - No interpretability libraries\n")
        return False


def test_file_structure():
    """Test that project structure is complete."""
    print("=" * 80)
    print("TEST 7: Project Structure")
    print("=" * 80)
    
    required_files = [
        "README.md",
        "requirements.txt",
        "config.py",
        "src/preprocessing.py",
        "src/models.py",
        "src/interpretability.py",
        "app/backend/main.py",
        "outputs/dataset_final.csv",
        "outputs/models/best_model_xgboost.pkl",
        "outputs/preprocessors/preprocessor.pkl",
        "docs/project_proposal.md",
        "docs/MODEL_CARD.md",
        "docs/API_DOCS.md",
        "docs/presentation/PRESENTATION_SLIDES.md",
        "tests/test_preprocessing.py",
        "tests/test_models.py",
        "tests/test_api.py",
        "Dockerfile.backend",
        "docker-compose.yml"
    ]
    
    missing = []
    present = []
    
    for file_path in required_files:
        if os.path.exists(file_path):
            present.append(file_path)
        else:
            missing.append(file_path)
    
    print(f"Files present: {len(present)}/{len(required_files)}")
    
    if missing:
        print(f"\n[WARN] Missing files ({len(missing)}):")
        for f in missing[:5]:  # Show first 5
            print(f"  - {f}")
        if len(missing) > 5:
            print(f"  ... and {len(missing) - 5} more")
    
    if len(present) >= len(required_files) * 0.8:  # 80% threshold
        print("[PASS] Test 7 PASSED (most files present)\n")
        return True
    else:
        print("[FAIL] Test 7 FAILED (too many missing files)\n")
        return False


def main():
    """Run all tests."""
    print("\n")
    print("*" + "=" * 78 + "*")
    print("|" + " " * 20 + "SYSTEM INTEGRATION TEST" + " " * 35 + "|")
    print("|" + " " * 15 + "Honey Bee Toxicity Prediction System" + " " * 27 + "|")
    print("*" + "=" * 78 + "*")
    print("\n")
    
    tests = [
        ("Model & Preprocessor Files", test_model_files),
        ("Dataset Loading", test_data_loading),
        ("Prediction Pipeline", test_prediction_pipeline),
        ("API Dependencies", test_api_imports),
        ("Visualization Libraries", test_visualization_libraries),
        ("Interpretability Libraries", test_interpretability_libraries),
        ("Project Structure", test_file_structure),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"[FAIL] Test '{test_name}' crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n")
    print("=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "[PASS]" if result else "[FAIL]"
        print(f"{status:8} | {test_name}")
    
    print("-" * 80)
    print(f"Total: {passed}/{total} tests passed ({passed/total*100:.1f}%)")
    
    if passed == total:
        print("\n[SUCCESS] ALL TESTS PASSED! System is ready for deployment.")
        return 0
    elif passed >= total * 0.7:
        print("\n[WARN] Most tests passed. System is mostly functional.")
        return 1
    else:
        print("\n[FAIL] Multiple tests failed. Review errors above.")
        return 2


if __name__ == "__main__":
    sys.exit(main())

