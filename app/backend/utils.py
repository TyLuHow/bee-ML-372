#!/usr/bin/env python3
"""
Backend Utility Functions
=========================

Shared utility functions for model loading, data validation,
and error handling across the backend.

Author: IME 372 Project Team
Date: November 2025
"""

import joblib
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional, List
import os
import json
from pathlib import Path


def load_model(path: str, model_name: str = "model") -> Any:
    """
    Load a joblib model with error handling.

    Args:
        path: Path to the model file
        model_name: Descriptive name for logging

    Returns:
        Loaded model object

    Raises:
        FileNotFoundError: If model file doesn't exist
        Exception: If model loading fails
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{model_name} not found at {path}")

    try:
        model = joblib.load(path)
        print(f"✓ {model_name} loaded from {path}")
        return model
    except Exception as e:
        raise Exception(f"Error loading {model_name} from {path}: {str(e)}")


def load_data(path: str, required_columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Load CSV data with validation.

    Args:
        path: Path to CSV file
        required_columns: List of required column names

    Returns:
        Loaded DataFrame

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If required columns are missing
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found at {path}")

    try:
        df = pd.read_csv(path)
        print(f"✓ Data loaded from {path} ({len(df)} rows)")

        # Validate required columns
        if required_columns:
            missing_cols = set(required_columns) - set(df.columns)
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")

        return df

    except Exception as e:
        raise Exception(f"Error loading data from {path}: {str(e)}")


def load_json(path: str, description: str = "JSON data") -> Dict[str, Any]:
    """
    Load JSON file with error handling.

    Args:
        path: Path to JSON file
        description: Description for logging

    Returns:
        Loaded dictionary

    Raises:
        FileNotFoundError: If file doesn't exist
        json.JSONDecodeError: If JSON is invalid
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{description} not found at {path}")

    try:
        with open(path, 'r') as f:
            data = json.load(f)
        print(f"✓ {description} loaded from {path}")
        return data
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(f"Invalid JSON in {path}: {str(e)}", e.doc, e.pos)


def save_json(data: Dict[str, Any], path: str, description: str = "data") -> None:
    """
    Save dictionary to JSON file with error handling.

    Args:
        data: Dictionary to save
        path: Path to save to
        description: Description for logging
    """
    try:
        # Create directory if needed
        os.makedirs(os.path.dirname(path), exist_ok=True)

        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"✓ {description} saved to {path}")
    except Exception as e:
        print(f"Warning: Could not save {description} to {path}: {e}")


def format_error(exception: Exception, context: str = "") -> Dict[str, str]:
    """
    Standardize error message formatting.

    Args:
        exception: Exception object
        context: Additional context about where error occurred

    Returns:
        Dictionary with formatted error message
    """
    error_msg = str(exception)
    if context:
        error_msg = f"{context}: {error_msg}"

    return {
        "error": error_msg,
        "error_type": type(exception).__name__
    }


def validate_input_features(
    data: pd.DataFrame,
    expected_features: List[str],
    add_missing: bool = True
) -> pd.DataFrame:
    """
    Validate and align input features with expected feature set.

    Args:
        data: Input DataFrame
        expected_features: List of expected feature names
        add_missing: If True, add missing columns with 0s

    Returns:
        DataFrame with aligned features

    Raises:
        ValueError: If required features are missing and add_missing=False
    """
    missing_features = set(expected_features) - set(data.columns)

    if missing_features:
        if add_missing:
            for col in missing_features:
                data[col] = 0
        else:
            raise ValueError(f"Missing required features: {missing_features}")

    # Reorder to match expected order
    return data[expected_features]


def calculate_confidence_interval(
    proportion: float,
    n: int,
    confidence: float = 0.95
) -> tuple:
    """
    Calculate Wilson score confidence interval for a proportion.

    Args:
        proportion: Observed proportion (0-1)
        n: Sample size
        confidence: Confidence level (default 0.95 for 95% CI)

    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    if n == 0:
        return (0.0, 0.0)

    # Z-score for confidence level
    z_scores = {0.90: 1.645, 0.95: 1.96, 0.99: 2.576}
    z = z_scores.get(confidence, 1.96)

    p = proportion
    denominator = 1 + z**2/n
    center = (p + z**2/(2*n)) / denominator
    margin = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denominator

    ci_lower = max(0, center - margin)
    ci_upper = min(1, center + margin)

    return (ci_lower, ci_upper)


def ensure_directory(path: str) -> None:
    """
    Ensure directory exists, create if it doesn't.

    Args:
        path: Directory path
    """
    Path(path).mkdir(parents=True, exist_ok=True)


def get_file_age(path: str) -> float:
    """
    Get age of file in seconds.

    Args:
        path: File path

    Returns:
        Age in seconds, or float('inf') if file doesn't exist
    """
    if not os.path.exists(path):
        return float('inf')

    import time
    return time.time() - os.path.getmtime(path)


def is_valid_smiles(smiles: str) -> bool:
    """
    Check if SMILES string is valid using RDKit.

    Args:
        smiles: SMILES string

    Returns:
        True if valid, False otherwise
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


def truncate_history(history: List[Dict], max_length: int = 100) -> List[Dict]:
    """
    Truncate history list to maximum length, keeping most recent.

    Args:
        history: List of history entries
        max_length: Maximum length to keep

    Returns:
        Truncated history list
    """
    if len(history) > max_length:
        return history[-max_length:]
    return history
