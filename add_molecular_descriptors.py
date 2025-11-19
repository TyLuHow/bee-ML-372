#!/usr/bin/env python3
"""
Add Molecular Descriptors to Dataset
=====================================

Reads the dataset, calculates RDKit molecular descriptors from SMILES,
and saves an enriched dataset with all features.

Author: IME 372 Project Team
Date: November 2025
"""

import pandas as pd
import os
from src.molecular_features import MolecularFeaturizer


def main():
    print("="*80)
    print("ADDING MOLECULAR DESCRIPTORS TO DATASET")
    print("="*80)
    
    # Load original dataset
    input_path = 'outputs/dataset_final.csv'
    output_path = 'outputs/dataset_with_descriptors.csv'
    
    print(f"\n1. Loading dataset from {input_path}...")
    df = pd.read_csv(input_path)
    print(f"   Loaded {len(df)} compounds with {len(df.columns)} columns")
    print(f"   Original columns: {df.columns.tolist()}")
    
    # Check for SMILES column
    if 'SMILES' not in df.columns:
        print("\nERROR: SMILES column not found in dataset!")
        return 1
    
    print(f"\n2. Calculating molecular descriptors from SMILES...")
    featurizer = MolecularFeaturizer()
    
    # Process all SMILES
    smiles_list = df['SMILES'].tolist()
    descriptors_df = featurizer.batch_smiles_to_dataframe(smiles_list)
    
    print(f"   Calculated descriptors for {len(descriptors_df)} compounds")
    print(f"   Descriptor columns: {descriptors_df.columns.tolist()}")
    
    # Merge descriptors with original data
    print(f"\n3. Merging descriptors with original data...")
    
    # Reset indices to align
    df_reset = df.reset_index(drop=True)
    descriptors_reset = descriptors_df.reset_index(drop=True)
    
    # Drop duplicate columns from descriptors (keep SMILES from original, drop from descriptors)
    cols_to_drop = [col for col in descriptors_reset.columns if col in df_reset.columns and col != 'index']
    if cols_to_drop:
        descriptors_reset = descriptors_reset.drop(columns=cols_to_drop)
    
    # Also drop the 'index' column from descriptors if it exists
    if 'index' in descriptors_reset.columns:
        descriptors_reset = descriptors_reset.drop(columns=['index'])
    
    # Combine
    df_enriched = pd.concat([df_reset, descriptors_reset], axis=1)
    
    print(f"   Enriched dataset shape: {df_enriched.shape}")
    print(f"   Total columns: {len(df_enriched.columns)}")
    
    # Save enriched dataset
    print(f"\n4. Saving enriched dataset to {output_path}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_enriched.to_csv(output_path, index=False)
    
    print(f"   OK Saved successfully!")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Original columns:  {len(df.columns)}")
    print(f"Descriptors added: {len(descriptors_df.columns)}")
    print(f"Final columns:     {len(df_enriched.columns)}")
    print(f"Total compounds:   {len(df_enriched)}")
    
    # Show new molecular descriptor columns
    descriptor_cols = [col for col in df_enriched.columns if col not in df.columns]
    print(f"\nNew molecular descriptor columns ({len(descriptor_cols)}):")
    for col in descriptor_cols[:10]:  # Show first 10
        print(f"  - {col}")
    if len(descriptor_cols) > 10:
        print(f"  ... and {len(descriptor_cols) - 10} more")
    
    print("\n" + "="*80)
    print("OK MOLECULAR DESCRIPTORS ADDED SUCCESSFULLY")
    print("="*80)
    print(f"\nOutput: {output_path}")
    print("You can now run:")
    print("  - python src/models.py")
    print("  - python run_all_analyses.py")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

