#!/usr/bin/env python3
"""
Test Scaffold-Based Splitting
==============================

Compare scaffold-based splitting with random splitting to assess
model generalization to novel chemical scaffolds.

Author: IME 372 Project Team
Date: November 2025
"""

from src.preprocessing import DataPreprocessor
from src.models import ModelTrainer
import pandas as pd
import warnings
warnings.filterwarnings('ignore')


def main():
    print("="*80)
    print("SCAFFOLD SPLITTING COMPARISON")
    print("="*80)
    
    # Load data
    df = pd.read_csv('outputs/dataset_with_descriptors.csv')
    print(f"\n[FILE] Loaded {len(df)} compounds")
    
    # Check if SMILES column exists
    if 'SMILES' not in df.columns:
        print("\nWARNING:  SMILES column not found in dataset")
        print("  Scaffold splitting requires SMILES strings")
        print("  Falling back to standard split comparison")
        return
    
    preprocessor = DataPreprocessor(random_state=42)
    
    # Prepare features
    feature_cols = [col for col in df.columns if col not in [
        'CID', 'CAS', 'name', 'SMILES', 'label', 'ppdb_level',
        'source', 'year', 'herbicide', 'fungicide', 'insecticide',
        'other_agrochemical', 'toxicity_type'
    ]]
    
    X = df[feature_cols + ['SMILES']].copy()
    y = df['label'].copy()
    
    print(f"Using {len(feature_cols)} features")
    
    # 1. Random split
    print("\n" + "="*80)
    print("1    RANDOM SPLIT")
    print("="*80)
    
    X_random = X.drop('SMILES', axis=1, errors='ignore')
    X_train_rand, X_val_rand, X_test_rand, y_train_rand, y_val_rand, y_test_rand = \
        preprocessor.split_data(X_random, y, test_size=0.2, val_size=0.1)
    
    # 2. Scaffold split
    print("\n" + "="*80)
    print("2    SCAFFOLD SPLIT")
    print("="*80)
    
    try:
        scaffold_splits = preprocessor.scaffold_split(X, y, smiles_col='SMILES')
        
        # Remove SMILES before training
        X_train_scaf = scaffold_splits['X_train'].drop('SMILES', axis=1, errors='ignore')
        X_val_scaf = scaffold_splits['X_val'].drop('SMILES', axis=1, errors='ignore')
        X_test_scaf = scaffold_splits['X_test'].drop('SMILES', axis=1, errors='ignore')
        y_train_scaf = scaffold_splits['y_train']
        y_val_scaf = scaffold_splits['y_val']
        y_test_scaf = scaffold_splits['y_test']
        
        # Train model on both splits
        trainer = ModelTrainer(random_state=42)
        
        print("\n" + "="*80)
        print("3    TRAINING ON RANDOM SPLIT")
        print("="*80)
        
        result_random = trainer.train_model(
            'xgboost',
            X_train_rand, y_train_rand,
            X_val_rand, y_val_rand,
            tune_hyperparams=False
        )
        
        print("\n" + "="*80)
        print("4    TRAINING ON SCAFFOLD SPLIT")
        print("="*80)
        
        result_scaffold = trainer.train_model(
            'xgboost',
            X_train_scaf, y_train_scaf,
            X_val_scaf, y_val_scaf,
            tune_hyperparams=False
        )
        
        # Compare results
        print("\n" + "="*80)
        print("RESULTS COMPARISON")
        print("="*80)
        
        print("\n[DATA] Validation Performance:")
        print(f"\n  Random Split:")
        print(f"    - Accuracy: {result_random['val_metrics']['accuracy']:.4f}")
        print(f"    - F1 Score: {result_random['val_metrics']['f1']:.4f}")
        print(f"    - ROC-AUC:  {result_random['val_metrics']['roc_auc']:.4f}")
        
        print(f"\n  Scaffold Split:")
        print(f"    - Accuracy: {result_scaffold['val_metrics']['accuracy']:.4f}")
        print(f"    - F1 Score: {result_scaffold['val_metrics']['f1']:.4f}")
        print(f"    - ROC-AUC:  {result_scaffold['val_metrics']['roc_auc']:.4f}")
        
        performance_drop = (
            (result_random['val_metrics']['f1'] - result_scaffold['val_metrics']['f1']) /
            result_random['val_metrics']['f1'] * 100
        )
        
        print(f"\n    Performance Drop (F1): {performance_drop:.1f}%")
        
        if performance_drop < 5:
            print(f"  [OK] Excellent generalization to novel scaffolds!")
        elif performance_drop < 10:
            print(f"  OK Good generalization to novel scaffolds")
        elif performance_drop < 20:
            print(f"  WARNING:  Moderate scaffold dependence detected")
        else:
            print(f"  [FAIL] Significant scaffold dependence - model may not generalize well")
        
        # Save results
        import json
        results = {
            'random_split': {
                'accuracy': float(result_random['val_metrics']['accuracy']),
                'f1': float(result_random['val_metrics']['f1']),
                'roc_auc': float(result_random['val_metrics']['roc_auc'])
            },
            'scaffold_split': {
                'accuracy': float(result_scaffold['val_metrics']['accuracy']),
                'f1': float(result_scaffold['val_metrics']['f1']),
                'roc_auc': float(result_scaffold['val_metrics']['roc_auc']),
                'n_scaffolds_train': scaffold_splits['n_scaffolds_train'],
                'n_scaffolds_val': scaffold_splits['n_scaffolds_val'],
                'n_scaffolds_test': scaffold_splits['n_scaffolds_test']
            },
            'performance_drop_pct': float(performance_drop)
        }
        
        import os
        os.makedirs('outputs/analysis', exist_ok=True)
        with open('outputs/analysis/scaffold_split_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n  OK Results saved to outputs/analysis/scaffold_split_results.json")
        
    except Exception as e:
        print(f"\nX Scaffold split failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("[OK] SCAFFOLD SPLIT TESTING COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()


