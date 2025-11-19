#!/usr/bin/env python3
"""
Data Preprocessing and Feature Engineering Module
==================================================

This module handles all data preprocessing steps including:
- Feature encoding (categorical variables)
- Feature scaling/normalization
- Feature selection
- Train/test splitting with stratification
- Handling class imbalance

Author: IME 372 Project Team
Date: November 2025
"""

import pandas as pd
import numpy as np
from typing import Tuple, List, Optional, Dict
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from imblearn.combine import SMOTETomek
import joblib
import warnings
warnings.filterwarnings('ignore')


class DataPreprocessor:
    """
    Comprehensive data preprocessing pipeline for ML model training.
    
    Features:
    - Categorical encoding
    - Feature scaling
    - Feature selection
    - Class imbalance handling
    - Reproducible preprocessing
    """
    
    def __init__(self, random_state: int = 42):
        """
        Initialize the preprocessor.
        
        Args:
            random_state: Random seed for reproducibility
        """
        self.random_state = random_state
        self.scaler = StandardScaler()
        self.feature_names = None
        self.label_encoders = {}
        self.onehot_encoder = None
        self.feature_selector = None
        self.selected_features = None
        
    def load_data(self, filepath: str) -> pd.DataFrame:
        """
        Load dataset from CSV file.
        
        Args:
            filepath: Path to CSV file
            
        Returns:
            Loaded dataframe
        """
        df = pd.read_csv(filepath)
        print(f"  Data loaded: {df.shape}")
        return df
    
    def prepare_features(
        self, 
        df: pd.DataFrame,
        target_col: str = 'label',
        exclude_cols: Optional[List[str]] = None
    ) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Prepare features by separating target and removing non-feature columns.
        
        Args:
            df: Input dataframe
            target_col: Name of target column
            exclude_cols: List of columns to exclude from features
            
        Returns:
            Tuple of (features dataframe, target series)
        """
        if exclude_cols is None:
            # Default columns to exclude (non-predictive features)
            exclude_cols = ['name', 'CID', 'CAS', 'SMILES', 'ppdb_level']
        
        # Separate target
        y = df[target_col].copy()
        
        # Remove target and excluded columns
        X = df.drop([target_col] + exclude_cols, axis=1, errors='ignore')
        
        print(f"  Features prepared: {X.shape}")
        print(f"  - Target: {target_col}")
        print(f"  - Features: {list(X.columns)}")
        
        return X, y
    
    def encode_categorical_features(
        self, 
        X: pd.DataFrame,
        categorical_cols: Optional[List[str]] = None,
        method: str = 'onehot'
    ) -> pd.DataFrame:
        """
        Encode categorical features.
        
        Args:
            X: Input features
            categorical_cols: List of categorical columns (auto-detected if None)
            method: Encoding method ('onehot' or 'label')
            
        Returns:
            Encoded features
        """
        X = X.copy()
        
        # Auto-detect categorical columns if not specified
        if categorical_cols is None:
            categorical_cols = X.select_dtypes(include=['object']).columns.tolist()
        
        if not categorical_cols:
            print("  No categorical features to encode")
            return X
        
        print(f"\nEncoding categorical features: {categorical_cols}")
        
        if method == 'onehot':
            # One-hot encoding
            X_encoded = pd.get_dummies(X, columns=categorical_cols, drop_first=True)
            print(f"  One-hot encoding complete: {X.shape}   {X_encoded.shape}")
            return X_encoded
        
        elif method == 'label':
            # Label encoding
            for col in categorical_cols:
                if col not in self.label_encoders:
                    self.label_encoders[col] = LabelEncoder()
                    X[col] = self.label_encoders[col].fit_transform(X[col])
                else:
                    X[col] = self.label_encoders[col].transform(X[col])
            print(f"  Label encoding complete")
            return X
        
        else:
            raise ValueError(f"Unknown encoding method: {method}")
    
    def scale_features(
        self, 
        X: pd.DataFrame, 
        fit: bool = True
    ) -> pd.DataFrame:
        """
        Scale numerical features using StandardScaler.
        
        Args:
            X: Input features
            fit: Whether to fit the scaler (True for training, False for test)
            
        Returns:
            Scaled features
        """
        if fit:
            X_scaled = self.scaler.fit_transform(X)
            print(f"  Features scaled (fitted)")
        else:
            X_scaled = self.scaler.transform(X)
            print(f"  Features scaled (transformed)")
        
        # Convert back to dataframe
        X_scaled = pd.DataFrame(X_scaled, columns=X.columns, index=X.index)
        
        return X_scaled
    
    def select_features(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        method: str = 'mutual_info',
        k: int = 20,
        fit: bool = True
    ) -> pd.DataFrame:
        """
        Select top k features based on statistical tests.
        
        Args:
            X: Input features
            y: Target variable
            method: Selection method ('f_classif' or 'mutual_info')
            k: Number of features to select (or 'all')
            fit: Whether to fit the selector
            
        Returns:
            Selected features
        """
        if k == 'all' or k >= X.shape[1]:
            print("  Keeping all features")
            return X
        
        if fit:
            # Choose scoring function
            if method == 'f_classif':
                score_func = f_classif
            elif method == 'mutual_info':
                score_func = mutual_info_classif
            else:
                raise ValueError(f"Unknown method: {method}")
            
            # Fit selector
            self.feature_selector = SelectKBest(score_func=score_func, k=k)
            X_selected = self.feature_selector.fit_transform(X, y)
            
            # Get selected feature names
            selected_mask = self.feature_selector.get_support()
            self.selected_features = X.columns[selected_mask].tolist()
            
            print(f"  Feature selection ({method}): {X.shape[1]}   {k} features")
            print(f"  Top features: {self.selected_features[:10]}")
            
        else:
            # Transform using fitted selector
            if self.feature_selector is None:
                raise ValueError("Feature selector not fitted!")
            X_selected = self.feature_selector.transform(X)
        
        # Convert back to dataframe
        X_selected = pd.DataFrame(
            X_selected, 
            columns=self.selected_features if fit else X.columns[:k],
            index=X.index
        )
        
        return X_selected
    
    def handle_imbalance(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        method: str = 'smote',
        sampling_strategy: str = 'auto'
    ) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Handle class imbalance using resampling techniques.
        
        Args:
            X: Input features
            y: Target variable
            method: Resampling method ('smote', 'undersample', 'smote_tomek', 'none')
            sampling_strategy: Sampling strategy ('auto', 'minority', 'not minority', 'all', or float)
            
        Returns:
            Tuple of (resampled features, resampled target)
        """
        if method == 'none':
            print("  No resampling applied")
            return X, y
        
        print(f"\nHandling class imbalance using: {method}")
        print(f"  Before: {y.value_counts().to_dict()}")
        
        if method == 'smote':
            sampler = SMOTE(random_state=self.random_state, sampling_strategy=sampling_strategy)
        elif method == 'undersample':
            sampler = RandomUnderSampler(random_state=self.random_state, sampling_strategy=sampling_strategy)
        elif method == 'smote_tomek':
            sampler = SMOTETomek(random_state=self.random_state, sampling_strategy=sampling_strategy)
        else:
            raise ValueError(f"Unknown resampling method: {method}")
        
        X_resampled, y_resampled = sampler.fit_resample(X, y)
        
        # Convert back to dataframe/series
        X_resampled = pd.DataFrame(X_resampled, columns=X.columns)
        y_resampled = pd.Series(y_resampled, name=y.name)
        
        print(f"  After:  {y_resampled.value_counts().to_dict()}")
        print(f"  Resampling complete: {len(y)}   {len(y_resampled)} samples")
        
        return X_resampled, y_resampled
    
    def split_data(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        test_size: float = 0.2,
        val_size: float = 0.1,
        stratify: bool = True
    ) -> Tuple:
        """
        Split data into train, validation, and test sets with stratification.
        
        Args:
            X: Features
            y: Target
            test_size: Proportion of data for test set
            val_size: Proportion of training data for validation set
            stratify: Whether to stratify splits
            
        Returns:
            Tuple of (X_train, X_val, X_test, y_train, y_val, y_test)
        """
        # First split: train+val vs test
        stratify_arg = y if stratify else None
        
        X_temp, X_test, y_temp, y_test = train_test_split(
            X, y, 
            test_size=test_size, 
            random_state=self.random_state,
            stratify=stratify_arg
        )
        
        # Second split: train vs val
        val_size_adjusted = val_size / (1 - test_size)
        stratify_arg_val = y_temp if stratify else None
        
        X_train, X_val, y_train, y_val = train_test_split(
            X_temp, y_temp,
            test_size=val_size_adjusted,
            random_state=self.random_state,
            stratify=stratify_arg_val
        )
        
        print(f"\n  Data split complete:")
        print(f"  Train: {len(X_train)} samples ({len(X_train)/len(X)*100:.1f}%)")
        print(f"  Val:   {len(X_val)} samples ({len(X_val)/len(X)*100:.1f}%)")
        print(f"  Test:  {len(X_test)} samples ({len(X_test)/len(X)*100:.1f}%)")
        
        if stratify:
            print(f"\n  Train class distribution: {y_train.value_counts().to_dict()}")
            print(f"  Val class distribution:   {y_val.value_counts().to_dict()}")
            print(f"  Test class distribution:  {y_test.value_counts().to_dict()}")
        
        return X_train, X_val, X_test, y_train, y_val, y_test
    
    def scaffold_split(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        smiles_col: str = 'SMILES',
        test_size: float = 0.2,
        val_size: float = 0.1
    ) -> Dict:
        """
        Split data by molecular scaffolds for structural diversity testing.
        
        Ensures train/val/test sets contain different core structures to test
        generalization to novel chemical scaffolds.
        
        Args:
            X: Feature dataframe (must include SMILES column)
            y: Target series
            smiles_col: Name of SMILES column
            test_size: Fraction for test set
            val_size: Fraction for validation set
            
        Returns:
            Dictionary with train/val/test splits and scaffold statistics
        """
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from collections import defaultdict
        
        print("\n  Performing scaffold-based split...")
        
        # Ensure SMILES column exists
        if smiles_col not in X.columns:
            raise ValueError(f"SMILES column '{smiles_col}' not found in dataframe")
        
        # Group molecules by scaffold
        scaffold_to_indices = defaultdict(list)
        invalid_smiles = []
        
        for idx, smiles in enumerate(X[smiles_col]):
            try:
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is None:
                    scaffold = "INVALID"
                    invalid_smiles.append(idx)
                else:
                    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
            except Exception as e:
                scaffold = "INVALID"
                invalid_smiles.append(idx)
            
            scaffold_to_indices[scaffold].append(idx)
        
        print(f"  - Unique scaffolds found: {len(scaffold_to_indices)}")
        print(f"  - Invalid SMILES: {len(invalid_smiles)}")
        
        # Sort scaffolds by size (largest first) for better distribution
        scaffolds = sorted(
            scaffold_to_indices.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )
        
        # Calculate split sizes
        n_total = len(X)
        n_test = int(n_total * test_size)
        n_val = int(n_total * val_size)
        n_train = n_total - n_test - n_val
        
        # Assign scaffolds to splits
        train_idx, val_idx, test_idx = [], [], []
        
        for scaffold, indices in scaffolds:
            if len(train_idx) < n_train:
                train_idx.extend(indices)
            elif len(val_idx) < n_val:
                val_idx.extend(indices)
            else:
                test_idx.extend(indices)
        
        # Count unique scaffolds per split
        train_scaffolds = set()
        val_scaffolds = set()
        test_scaffolds = set()
        
        for scaffold, indices in scaffolds:
            if any(idx in train_idx for idx in indices):
                train_scaffolds.add(scaffold)
            if any(idx in val_idx for idx in indices):
                val_scaffolds.add(scaffold)
            if any(idx in test_idx for idx in indices):
                test_scaffolds.add(scaffold)
        
        # Create splits
        splits = {
            'X_train': X.iloc[train_idx].reset_index(drop=True),
            'y_train': y.iloc[train_idx].reset_index(drop=True),
            'X_val': X.iloc[val_idx].reset_index(drop=True),
            'y_val': y.iloc[val_idx].reset_index(drop=True),
            'X_test': X.iloc[test_idx].reset_index(drop=True),
            'y_test': y.iloc[test_idx].reset_index(drop=True),
            'split_type': 'scaffold',
            'n_scaffolds_train': len(train_scaffolds),
            'n_scaffolds_val': len(val_scaffolds),
            'n_scaffolds_test': len(test_scaffolds),
            'scaffold_overlap': {
                'train_val': len(train_scaffolds & val_scaffolds),
                'train_test': len(train_scaffolds & test_scaffolds),
                'val_test': len(val_scaffolds & test_scaffolds)
            }
        }
        
        print(f"\n    Scaffold split created:")
        print(f"    Train: {len(train_idx)} samples, {len(train_scaffolds)} scaffolds")
        print(f"    Val:   {len(val_idx)} samples, {len(val_scaffolds)} scaffolds")
        print(f"    Test:  {len(test_idx)} samples, {len(test_scaffolds)} scaffolds")
        print(f"\n  Scaffold overlap:")
        print(f"    Train-Val:  {splits['scaffold_overlap']['train_val']}")
        print(f"    Train-Test: {splits['scaffold_overlap']['train_test']}")
        print(f"    Val-Test:   {splits['scaffold_overlap']['val_test']}")
        
        return splits
    
    def save_preprocessor(self, filepath: str):
        """
        Save preprocessor instance for later use.
        
        Args:
            filepath: Path to save preprocessor
        """
        import os
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        # Save the entire DataPreprocessor instance (not dict)
        joblib.dump(self, filepath)
        print(f"  Preprocessor saved: {filepath}")
    
    @classmethod
    def load_preprocessor(cls, filepath: str) -> 'DataPreprocessor':
        """
        Load saved preprocessor instance.
        
        Args:
            filepath: Path to saved preprocessor
            
        Returns:
            Loaded DataPreprocessor instance
        """
        preprocessor = joblib.load(filepath)
        
        # Verify it's a DataPreprocessor instance
        if not isinstance(preprocessor, cls):
            raise TypeError(f"Loaded object is not a DataPreprocessor instance, got {type(preprocessor)}")
        
        print(f"  Preprocessor loaded: {filepath}")
        return preprocessor


def create_preprocessing_pipeline(
    data_path: str,
    target_col: str = 'label',
    test_size: float = 0.2,
    val_size: float = 0.1,
    apply_scaling: bool = True,
    apply_feature_selection: bool = False,
    n_features: int = 20,
    handle_imbalance: bool = True,
    imbalance_method: str = 'smote',
    save_path: str = 'outputs/preprocessors/preprocessor.pkl'
) -> Tuple:
    """
    Complete preprocessing pipeline from raw data to model-ready splits.
    
    Args:
        data_path: Path to dataset CSV
        target_col: Target column name
        test_size: Test set proportion
        val_size: Validation set proportion
        apply_scaling: Whether to scale features
        apply_feature_selection: Whether to apply feature selection
        n_features: Number of features to select (if feature selection applied)
        handle_imbalance: Whether to handle class imbalance
        imbalance_method: Method for handling imbalance ('smote', 'undersample', etc.)
        save_path: Path to save preprocessor
        
    Returns:
        Tuple of (X_train, X_val, X_test, y_train, y_val, y_test, preprocessor)
    """
    print("="*80)
    print("DATA PREPROCESSING PIPELINE")
    print("="*80)
    
    # Initialize preprocessor
    preprocessor = DataPreprocessor()
    
    # Load data
    df = preprocessor.load_data(data_path)
    
    # Prepare features
    X, y = preprocessor.prepare_features(df, target_col=target_col)
    
    # Encode categorical features
    X = preprocessor.encode_categorical_features(X, method='onehot')
    
    # Split data first (before scaling/resampling to avoid data leakage)
    X_train, X_val, X_test, y_train, y_val, y_test = preprocessor.split_data(
        X, y, test_size=test_size, val_size=val_size, stratify=True
    )
    
    # Scale features (fit on train only)
    if apply_scaling:
        X_train = preprocessor.scale_features(X_train, fit=True)
        X_val = preprocessor.scale_features(X_val, fit=False)
        X_test = preprocessor.scale_features(X_test, fit=False)
    
    # Feature selection (fit on train only)
    if apply_feature_selection:
        X_train = preprocessor.select_features(X_train, y_train, k=n_features, fit=True)
        X_val = preprocessor.select_features(X_val, y_val, k=n_features, fit=False)
        X_test = preprocessor.select_features(X_test, y_test, k=n_features, fit=False)
    
    # Handle class imbalance (only on training data)
    if handle_imbalance:
        X_train, y_train = preprocessor.handle_imbalance(
            X_train, y_train, method=imbalance_method
        )
    
    # Save preprocessor
    if save_path:
        preprocessor.save_preprocessor(save_path)
    
    print("\n" + "="*80)
    print("PREPROCESSING COMPLETE")
    print("="*80)
    print(f"Final shapes:")
    print(f"  X_train: {X_train.shape}")
    print(f"  X_val:   {X_val.shape}")
    print(f"  X_test:  {X_test.shape}")
    
    return X_train, X_val, X_test, y_train, y_val, y_test, preprocessor


if __name__ == "__main__":
    # Example usage
    X_train, X_val, X_test, y_train, y_val, y_test, preprocessor = create_preprocessing_pipeline(
        data_path='data/raw/dataset_with_descriptors.csv',
        apply_scaling=True,
        apply_feature_selection=False,
        handle_imbalance=True,
        imbalance_method='smote'
    )
    
    print("\n  Preprocessing pipeline executed successfully!")

