#!/usr/bin/env python3
"""
Model Training and Evaluation Module
=====================================

This module implements multiple ML algorithms for binary classification:
- Logistic Regression (baseline)
- Random Forest
- XGBoost
- LightGBM
- Support Vector Machine
- Multi-Layer Perceptron

Features:
- Hyperparameter tuning with GridSearchCV/RandomizedSearchCV
- Cross-validation with stratification
- Model comparison and selection
- Performance metrics tracking
- Model persistence

Author: IME 372 Project Team
Date: November 2025
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
import json
import time
from datetime import datetime
import joblib
import warnings
warnings.filterwarnings('ignore')

# Scikit-learn
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, StratifiedKFold, cross_val_score
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, confusion_matrix, classification_report,
    roc_curve, precision_recall_curve
)

# Gradient boosting libraries
import xgboost as xgb
import lightgbm as lgb

# Import preprocessing
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from preprocessing import create_preprocessing_pipeline


class ModelTrainer:
    """
    Unified model training and evaluation framework.
    """
    
    def __init__(self, random_state: int = 42, n_jobs: int = -1):
        """
        Initialize the model trainer.
        
        Args:
            random_state: Random seed for reproducibility
            n_jobs: Number of parallel jobs (-1 uses all cores)
        """
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.models = {}
        self.best_model = None
        self.results = {}
        
    def get_model(self, model_name: str, **kwargs) -> Any:
        """
        Get a model instance by name with optional custom parameters.
        
        Args:
            model_name: Name of the model
            **kwargs: Additional model parameters
            
        Returns:
            Instantiated model
        """
        models_dict = {
            'logistic': LogisticRegression(
                random_state=self.random_state,
                max_iter=1000,
                n_jobs=self.n_jobs,
                **kwargs
            ),
            'random_forest': RandomForestClassifier(
                random_state=self.random_state,
                n_jobs=self.n_jobs,
                **kwargs
            ),
            'xgboost': xgb.XGBClassifier(
                random_state=self.random_state,
                n_jobs=self.n_jobs,
                eval_metric='logloss',
                **kwargs
            ),
            'lightgbm': lgb.LGBMClassifier(
                random_state=self.random_state,
                n_jobs=self.n_jobs,
                verbose=-1,
                **kwargs
            ),
            'svm': SVC(
                random_state=self.random_state,
                probability=True,
                **kwargs
            ),
            'mlp': MLPClassifier(
                random_state=self.random_state,
                max_iter=500,
                **kwargs
            )
        }
        
        if model_name not in models_dict:
            raise ValueError(f"Unknown model: {model_name}. Choose from {list(models_dict.keys())}")
        
        return models_dict[model_name]
    
    def get_param_grid(self, model_name: str, search_type: str = 'grid') -> Dict:
        """
        Get hyperparameter grid for a model.
        
        Args:
            model_name: Name of the model
            search_type: 'grid' for GridSearch or 'random' for RandomizedSearch
            
        Returns:
            Parameter grid dictionary
        """
        if search_type == 'grid':
            # Smaller grids for GridSearch (faster)
            param_grids = {
                'logistic': {
                    'C': [0.1, 1.0, 10.0],
                    'penalty': ['l2'],
                    'solver': ['lbfgs', 'liblinear']
                },
                'random_forest': {
                    'n_estimators': [100, 200],
                    'max_depth': [10, 20, None],
                    'min_samples_split': [2, 5],
                    'min_samples_leaf': [1, 2],
                    'class_weight': ['balanced', None]
                },
                'xgboost': {
                    'n_estimators': [100, 200],
                    'max_depth': [3, 5, 7],
                    'learning_rate': [0.01, 0.1],
                    'subsample': [0.8, 1.0],
                    'colsample_bytree': [0.8, 1.0],
                    'scale_pos_weight': [1, 2.5]
                },
                'lightgbm': {
                    'n_estimators': [100, 200],
                    'max_depth': [3, 5, 7],
                    'learning_rate': [0.01, 0.1],
                    'num_leaves': [31, 63],
                    'subsample': [0.8, 1.0],
                    'colsample_bytree': [0.8, 1.0],
                    'class_weight': ['balanced', None]
                },
                'svm': {
                    'C': [0.1, 1.0, 10.0],
                    'kernel': ['rbf', 'linear'],
                    'gamma': ['scale', 'auto'],
                    'class_weight': ['balanced', None]
                },
                'mlp': {
                    'hidden_layer_sizes': [(50,), (100,), (50, 50)],
                    'activation': ['relu', 'tanh'],
                    'alpha': [0.0001, 0.001],
                    'learning_rate': ['constant', 'adaptive']
                }
            }
        else:
            # Larger ranges for RandomizedSearch
            param_grids = {
                'logistic': {
                    'C': [0.01, 0.1, 1.0, 10.0, 100.0],
                    'penalty': ['l1', 'l2'],
                    'solver': ['lbfgs', 'liblinear', 'saga']
                },
                'random_forest': {
                    'n_estimators': [50, 100, 200, 300],
                    'max_depth': [5, 10, 20, 30, None],
                    'min_samples_split': [2, 5, 10],
                    'min_samples_leaf': [1, 2, 4],
                    'max_features': ['sqrt', 'log2'],
                    'class_weight': ['balanced', 'balanced_subsample', None]
                },
                'xgboost': {
                    'n_estimators': [50, 100, 200, 300],
                    'max_depth': [3, 5, 7, 10],
                    'learning_rate': [0.01, 0.05, 0.1, 0.2],
                    'subsample': [0.6, 0.8, 1.0],
                    'colsample_bytree': [0.6, 0.8, 1.0],
                    'gamma': [0, 0.1, 0.2],
                    'scale_pos_weight': [1, 2, 2.5, 3]
                },
                'lightgbm': {
                    'n_estimators': [50, 100, 200, 300],
                    'max_depth': [3, 5, 7, 10, -1],
                    'learning_rate': [0.01, 0.05, 0.1, 0.2],
                    'num_leaves': [31, 63, 127],
                    'subsample': [0.6, 0.8, 1.0],
                    'colsample_bytree': [0.6, 0.8, 1.0],
                    'reg_alpha': [0, 0.1, 0.5],
                    'reg_lambda': [0, 0.1, 0.5],
                    'class_weight': ['balanced', None]
                },
                'svm': {
                    'C': [0.1, 1.0, 10.0, 100.0],
                    'kernel': ['rbf', 'linear', 'poly'],
                    'gamma': ['scale', 'auto', 0.001, 0.01],
                    'class_weight': ['balanced', None]
                },
                'mlp': {
                    'hidden_layer_sizes': [(50,), (100,), (50, 50), (100, 50), (100, 100)],
                    'activation': ['relu', 'tanh'],
                    'alpha': [0.0001, 0.001, 0.01],
                    'learning_rate': ['constant', 'invscaling', 'adaptive'],
                    'batch_size': [32, 64, 128]
                }
            }
        
        return param_grids.get(model_name, {})
    
    def train_model(
        self,
        model_name: str,
        X_train: pd.DataFrame,
        y_train: pd.Series,
        X_val: Optional[pd.DataFrame] = None,
        y_val: Optional[pd.Series] = None,
        tune_hyperparams: bool = True,
        search_type: str = 'grid',
        cv_folds: int = 5,
        scoring: str = 'f1'
    ) -> Dict[str, Any]:
        """
        Train a single model with optional hyperparameter tuning.
        
        Args:
            model_name: Name of the model to train
            X_train: Training features
            y_train: Training target
            X_val: Validation features (optional)
            y_val: Validation target (optional)
            tune_hyperparams: Whether to perform hyperparameter tuning
            search_type: 'grid' or 'random' search
            cv_folds: Number of cross-validation folds
            scoring: Scoring metric for tuning
            
        Returns:
            Dictionary with trained model and metrics
        """
        print(f"\n{'='*80}")
        print(f"Training: {model_name.upper()}")
        print(f"{'='*80}")
        
        start_time = time.time()
        
        # Get base model
        model = self.get_model(model_name)
        
        # Hyperparameter tuning
        if tune_hyperparams:
            print(f"Performing {search_type} search for hyperparameter tuning...")
            param_grid = self.get_param_grid(model_name, search_type)
            
            if not param_grid:
                print(f"No param grid defined for {model_name}, using default parameters")
                best_model = model
            else:
                cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
                
                if search_type == 'grid':
                    search = GridSearchCV(
                        model,
                        param_grid,
                        cv=cv,
                        scoring=scoring,
                        n_jobs=self.n_jobs,
                        verbose=1
                    )
                else:  # random
                    search = RandomizedSearchCV(
                        model,
                        param_grid,
                        n_iter=20,
                        cv=cv,
                        scoring=scoring,
                        n_jobs=self.n_jobs,
                        random_state=self.random_state,
                        verbose=1
                    )
                
                search.fit(X_train, y_train)
                best_model = search.best_estimator_
                
                print(f"\nBest parameters: {search.best_params_}")
                print(f"Best CV score: {search.best_score_:.4f}")
        else:
            # Train with default parameters
            print("Training with default parameters...")
            best_model = model
            best_model.fit(X_train, y_train)
        
        training_time = time.time() - start_time
        
        # Evaluate on training set
        y_train_pred = best_model.predict(X_train)
        y_train_proba = best_model.predict_proba(X_train)[:, 1]
        
        train_metrics = self._calculate_metrics(y_train, y_train_pred, y_train_proba, "Training")
        
        # Evaluate on validation set if provided
        val_metrics = {}
        if X_val is not None and y_val is not None:
            y_val_pred = best_model.predict(X_val)
            y_val_proba = best_model.predict_proba(X_val)[:, 1]
            val_metrics = self._calculate_metrics(y_val, y_val_pred, y_val_proba, "Validation")
        
        # Store model and results
        result = {
            'model': best_model,
            'model_name': model_name,
            'train_metrics': train_metrics,
            'val_metrics': val_metrics,
            'training_time': training_time,
            'timestamp': datetime.now().isoformat()
        }
        
        self.models[model_name] = best_model
        self.results[model_name] = result
        
        print(f"\n  {model_name} training complete (Time: {training_time:.2f}s)")
        
        return result
    
    def _calculate_metrics(
        self,
        y_true: pd.Series,
        y_pred: np.ndarray,
        y_proba: np.ndarray,
        dataset_name: str = ""
    ) -> Dict[str, float]:
        """
        Calculate comprehensive classification metrics.
        
        Args:
            y_true: True labels
            y_pred: Predicted labels
            y_proba: Prediction probabilities
            dataset_name: Name of dataset (for display)
            
        Returns:
            Dictionary of metrics
        """
        metrics = {
            'accuracy': accuracy_score(y_true, y_pred),
            'precision': precision_score(y_true, y_pred, zero_division=0),
            'recall': recall_score(y_true, y_pred, zero_division=0),
            'f1': f1_score(y_true, y_pred, zero_division=0),
            'roc_auc': roc_auc_score(y_true, y_proba)
        }
        
        cm = confusion_matrix(y_true, y_pred)
        metrics['confusion_matrix'] = cm.tolist()
        
        if dataset_name:
            print(f"\n{dataset_name} Metrics:")
            print(f"  Accuracy:  {metrics['accuracy']:.4f}")
            print(f"  Precision: {metrics['precision']:.4f}")
            print(f"  Recall:    {metrics['recall']:.4f}")
            print(f"  F1 Score:  {metrics['f1']:.4f}")
            print(f"  ROC-AUC:   {metrics['roc_auc']:.4f}")
            print(f"\n  Confusion Matrix:")
            print(f"    [[TN={cm[0,0]}, FP={cm[0,1]}]")
            print(f"     [FN={cm[1,0]}, TP={cm[1,1]}]]")
        
        return metrics
    
    def train_all_models(
        self,
        X_train: pd.DataFrame,
        y_train: pd.Series,
        X_val: pd.DataFrame,
        y_val: pd.Series,
        models_to_train: Optional[List[str]] = None,
        tune_hyperparams: bool = True
    ) -> Dict[str, Dict]:
        """
        Train multiple models and compare performance.
        
        Args:
            X_train: Training features
            y_train: Training target
            X_val: Validation features
            y_val: Validation target
            models_to_train: List of model names (trains all if None)
            tune_hyperparams: Whether to tune hyperparameters
            
        Returns:
            Dictionary of all model results
        """
        if models_to_train is None:
            models_to_train = ['logistic', 'random_forest', 'xgboost', 'lightgbm', 'svm', 'mlp']
        
        print("="*80)
        print("TRAINING MULTIPLE MODELS")
        print("="*80)
        print(f"Models to train: {models_to_train}")
        print(f"Hyperparameter tuning: {tune_hyperparams}")
        
        all_results = {}
        
        for model_name in models_to_train:
            try:
                result = self.train_model(
                    model_name,
                    X_train,
                    y_train,
                    X_val,
                    y_val,
                    tune_hyperparams=tune_hyperparams,
                    search_type='grid'
                )
                all_results[model_name] = result
            except Exception as e:
                print(f"\n  Error training {model_name}: {str(e)}")
                continue
        
        # Create comparison table
        self._print_model_comparison(all_results)
        
        # Select best model
        self._select_best_model(all_results)
        
        return all_results
    
    def _print_model_comparison(self, results: Dict[str, Dict]):
        """Print a comparison table of all models."""
        print("\n" + "="*80)
        print("MODEL COMPARISON")
        print("="*80)
        
        comparison_data = []
        for model_name, result in results.items():
            val_metrics = result.get('val_metrics', {})
            if val_metrics:
                comparison_data.append({
                    'Model': model_name,
                    'Accuracy': f"{val_metrics['accuracy']:.4f}",
                    'Precision': f"{val_metrics['precision']:.4f}",
                    'Recall': f"{val_metrics['recall']:.4f}",
                    'F1 Score': f"{val_metrics['f1']:.4f}",
                    'ROC-AUC': f"{val_metrics['roc_auc']:.4f}",
                    'Time (s)': f"{result['training_time']:.2f}"
                })
        
        if comparison_data:
            df_comparison = pd.DataFrame(comparison_data)
            print(df_comparison.to_string(index=False))
        else:
            print("No validation metrics available for comparison")
    
    def _select_best_model(self, results: Dict[str, Dict], metric: str = 'f1'):
        """Select the best model based on validation metric."""
        best_score = -1
        best_model_name = None
        
        for model_name, result in results.items():
            val_metrics = result.get('val_metrics', {})
            if val_metrics and val_metrics.get(metric, 0) > best_score:
                best_score = val_metrics[metric]
                best_model_name = model_name
        
        if best_model_name:
            self.best_model = self.models[best_model_name]
            print(f"\n  Best model: {best_model_name} (Validation {metric.upper()}: {best_score:.4f})")
        else:
            print("\n  Could not select best model")
    
    def evaluate_on_test(
        self,
        X_test: pd.DataFrame,
        y_test: pd.Series,
        model_name: Optional[str] = None
    ) -> Dict[str, float]:
        """
        Evaluate model on test set.
        
        Args:
            X_test: Test features
            y_test: Test target
            model_name: Name of model to evaluate (uses best model if None)
            
        Returns:
            Test metrics
        """
        if model_name is None:
            if self.best_model is None:
                raise ValueError("No best model selected!")
            model = self.best_model
            model_name = "Best Model"
        else:
            if model_name not in self.models:
                raise ValueError(f"Model {model_name} not found!")
            model = self.models[model_name]
        
        print(f"\n{'='*80}")
        print(f"TEST SET EVALUATION: {model_name}")
        print(f"{'='*80}")
        
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:, 1]
        
        test_metrics = self._calculate_metrics(y_test, y_pred, y_proba, "Test")
        
        # Print classification report
        print(f"\nClassification Report:")
        print(classification_report(y_test, y_pred, target_names=['Non-toxic', 'Toxic']))
        
        return test_metrics
    
    def save_model(self, model_name: str, filepath: str):
        """Save a trained model."""
        if model_name not in self.models:
            raise ValueError(f"Model {model_name} not found!")
        
        joblib.dump(self.models[model_name], filepath)
        print(f"  Model saved: {filepath}")
    
    def save_results(self, filepath: str):
        """Save training results to JSON."""
        # Convert results to JSON-serializable format
        results_json = {}
        for model_name, result in self.results.items():
            results_json[model_name] = {
                'model_name': result['model_name'],
                'train_metrics': result['train_metrics'],
                'val_metrics': result['val_metrics'],
                'training_time': result['training_time'],
                'timestamp': result['timestamp']
            }
        
        with open(filepath, 'w') as f:
            json.dump(results_json, f, indent=2)
        
        print(f"  Results saved: {filepath}")


def main():
    """Main function to run the complete modeling pipeline."""
    print("="*80)
    print("PHASE 3: MODEL DEVELOPMENT & SELECTION")
    print("="*80)
    
    # Load preprocessed data
    print("\n1. Loading preprocessed data...")
    X_train, X_val, X_test, y_train, y_val, y_test, preprocessor = create_preprocessing_pipeline(
        data_path='data/raw/dataset_with_descriptors.csv',
        apply_scaling=True,
        apply_feature_selection=False,
        handle_imbalance=True,
        imbalance_method='smote'
    )
    
    # Initialize trainer
    trainer = ModelTrainer(random_state=42)
    
    # Train all models
    print("\n2. Training models...")
    all_results = trainer.train_all_models(
        X_train, y_train, X_val, y_val,
        models_to_train=['logistic', 'random_forest', 'xgboost', 'lightgbm'],
        tune_hyperparams=True
    )
    
    # Evaluate best model on test set
    print("\n3. Evaluating best model on test set...")
    test_metrics = trainer.evaluate_on_test(X_test, y_test)
    
    # Save best model
    print("\n4. Saving best model...")
    best_model_name = [name for name, result in all_results.items() 
                       if result['model'] == trainer.best_model][0]
    trainer.save_model(best_model_name, f'outputs/models/best_model_{best_model_name}.pkl')
    
    # Save results
    trainer.save_results('outputs/metrics/training_results.json')
    
    print("\n" + "="*80)
    print("  PHASE 3 COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()

