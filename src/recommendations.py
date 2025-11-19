#!/usr/bin/env python3
"""
Alternative Compound Recommendation System
===========================================

Finds safer alternative pesticides with similar molecular properties
using K-Nearest Neighbors in molecular descriptor space.

Author: IME 372 Project Team
Date: November 2025
"""

from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
from typing import List, Dict, Optional
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
import warnings
warnings.filterwarnings('ignore')


class CompoundRecommender:
    """
    Recommend safer alternative compounds with similar properties.
    
    Uses K-Nearest Neighbors in molecular descriptor space to find
    non-toxic alternatives to toxic compounds.
    """
    
    def __init__(
        self, 
        df: pd.DataFrame, 
        feature_cols: List[str],
        n_neighbors: int = 10
    ):
        """
        Initialize recommender with dataset.
        
        Args:
            df: DataFrame with compounds and toxicity labels
            feature_cols: List of molecular descriptor column names
            n_neighbors: Number of neighbors to consider
        """
        self.df = df.copy()
        self.feature_cols = feature_cols
        self.n_neighbors = n_neighbors
        
        # Separate toxic and non-toxic compounds
        self.toxic_df = df[df['label'] == 1].reset_index(drop=True)
        self.safe_df = df[df['label'] == 0].reset_index(drop=True)
        
        print(f"[DATA] Recommender initialized:")
        print(f"  - Total compounds: {len(df)}")
        print(f"  - Toxic compounds: {len(self.toxic_df)}")
        print(f"  - Safe compounds: {len(self.safe_df)}")
        print(f"  - Features: {len(feature_cols)}")
        
        # Scale features
        self.scaler = StandardScaler()
        safe_features = self.safe_df[feature_cols].fillna(0).values
        self.scaler.fit(safe_features)
        safe_features_scaled = self.scaler.transform(safe_features)
        
        # Fit KNN on safe compounds
        self.knn = NearestNeighbors(
            n_neighbors=min(n_neighbors, len(self.safe_df)),
            metric='euclidean',
            n_jobs=-1
        )
        self.knn.fit(safe_features_scaled)
        
        print(f"  OK KNN model trained on {len(self.safe_df)} safe compounds")
    
    def find_alternatives(
        self, 
        compound_id: str,
        id_col: str = 'CID',
        n_alternatives: int = 5
    ) -> Optional[pd.DataFrame]:
        """
        Find non-toxic alternatives similar to a compound.
        
        Args:
            compound_id: Identifier of compound to find alternatives for
            id_col: Column name for identifier (default: 'CID')
            n_alternatives: Number of alternatives to return
            
        Returns:
            DataFrame with alternatives ranked by similarity, or None if not found
        """
        # Find compound
        compound = self.df[self.df[id_col].astype(str) == str(compound_id)]
        
        if compound.empty:
            print(f"  X Compound {compound_id} not found")
            return None
        
        # Get features
        features = compound[self.feature_cols].fillna(0).values
        features_scaled = self.scaler.transform(features)
        
        # Find nearest safe neighbors
        n_request = min(n_alternatives, len(self.safe_df))
        distances, indices = self.knn.kneighbors(features_scaled, n_neighbors=n_request)
        
        # Get alternative compounds
        alternatives = self.safe_df.iloc[indices[0]].copy()
        alternatives['similarity_distance'] = distances[0]
        alternatives['similarity_score'] = 1 / (1 + distances[0])  # Convert to 0-1
        alternatives['rank'] = range(1, len(alternatives) + 1)
        
        # Add original compound info
        alternatives['original_name'] = compound['name'].values[0]
        alternatives['original_id'] = str(compound_id)
        alternatives['original_is_toxic'] = compound['label'].values[0] == 1
        
        return alternatives
    
    def batch_recommend(
        self, 
        save_path: str = 'outputs/analysis/alternatives.csv',
        n_alternatives: int = 3
    ) -> pd.DataFrame:
        """
        Generate recommendations for all toxic compounds.
        
        Args:
            save_path: Path to save recommendations
            n_alternatives: Number of alternatives per compound
            
        Returns:
            DataFrame with all recommendations
        """
        print("\n[PROCESS] Generating recommendations for all toxic compounds...")
        
        recommendations = []
        
        for idx, row in self.toxic_df.iterrows():
            cid = row.get('CID', row.get('name', f'compound_{idx}'))
            alts = self.find_alternatives(str(cid), n_alternatives=n_alternatives)
            
            if alts is not None:
                # Select key columns
                cols_to_keep = ['name', 'CID', 'rank', 'similarity_score', 'similarity_distance',
                               'original_name', 'original_id'] + self.feature_cols
                
                # Only keep columns that exist
                cols_available = [c for c in cols_to_keep if c in alts.columns]
                alt_subset = alts[cols_available].copy()
                
                recommendations.append(alt_subset)
        
        if recommendations:
            all_recommendations = pd.concat(recommendations, ignore_index=True)
            
            # Save
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            all_recommendations.to_csv(save_path, index=False)
            
            print(f"  OK Generated {len(recommendations)} recommendation sets")
            print(f"  OK Total alternatives: {len(all_recommendations)}")
            print(f"  OK Saved to {save_path}")
            
            return all_recommendations
        
        return pd.DataFrame()
    
    def create_summary_report(
        self, 
        recommendations_df: pd.DataFrame,
        output_path: str = 'outputs/ALTERNATIVES_REPORT.md'
    ):
        """
        Generate markdown report of recommendations.
        
        Args:
            recommendations_df: Output from batch_recommend()
            output_path: Path to save report
        """
        n_toxic = recommendations_df['original_id'].nunique()
        avg_similarity = recommendations_df.groupby('original_id')['similarity_score'].mean().mean()
        
        # Get top recommendations
        top_matches = recommendations_df[recommendations_df['rank'] == 1].nlargest(10, 'similarity_score')
        
        report = f"""# Alternative Compound Recommendations Report
**Generated**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Executive Summary

- **Toxic Compounds Analyzed**: {n_toxic}
- **Total Alternatives Generated**: {len(recommendations_df)}
- **Average Similarity Score**: {avg_similarity:.3f}
- **Recommendation Method**: K-Nearest Neighbors in molecular descriptor space

## Top 10 Best Alternatives (Rank 1, Highest Similarity)

| Original Compound | Alternative | Similarity Score | CID |
|-------------------|-------------|------------------|-----|
"""
        
        for idx, row in top_matches.iterrows():
            report += f"| {row['original_name']} | {row['name']} | {row['similarity_score']:.3f} | {row['CID']} |\n"
        
        report += """

## Methodology

### K-Nearest Neighbors Approach

The recommendation system uses K-Nearest Neighbors (KNN) to find structurally similar non-toxic compounds:

1. **Feature Space**: Molecular descriptors (MW, LogP, H-donors, etc.)
2. **Scaling**: StandardScaler normalization
3. **Distance Metric**: Euclidean distance
4. **Similarity Score**: 1 / (1 + distance)

### Advantages

- Fast retrieval of similar compounds
- Interpretable based on property similarity
- No training required (instance-based)

### Limitations

- Assumes toxicity is related to property similarity
- May miss compounds with different properties but similar safety profile
- Limited to compounds in training set

## Use Cases

### 1. Formulation Development

When developing new pesticide formulations, use this system to:
- Identify safer alternatives to toxic active ingredients
- Maintain similar physicochemical properties
- Reduce environmental impact

### 2. Risk Assessment

For regulatory review:
- Compare new compounds to known safe alternatives
- Assess if similar compounds have toxicity concerns
- Support read-across approaches

### 3. Green Chemistry

In sustainable pesticide design:
- Replace high-risk compounds with safer analogs
- Maintain efficacy through property matching
- Guide structure optimization

## API Usage

```python
# Get alternatives via API
GET /recommend/alternatives/CID12345?n_alternatives=5

# Response:
{
  "original_id": "CID12345",
  "original_name": "Compound Name",
  "n_alternatives": 5,
  "alternatives": [
    {
      "rank": 1,
      "name": "Alternative 1",
      "cid": "CID67890",
      "similarity_score": 0.92,
      "similarity_distance": 2.34
    },
    ...
  ]
}
```

---

*Report generated by ApisTox Alternative Compound Recommender*
"""
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(report)
        
        print(f"  OK Report saved to {output_path}")


# CLI interface
if __name__ == "__main__":
    print("="*80)
    print("ALTERNATIVE COMPOUND RECOMMENDER")
    print("="*80)
    
    # Load dataset
    df = pd.read_csv('outputs/dataset_with_descriptors.csv')
    print(f"\n[FILE] Loaded {len(df)} compounds")
    
    # Define features (molecular descriptors)
    feature_cols = [
        'MolecularWeight', 'LogP', 'NumHDonors', 'NumHAcceptors',
        'NumRotatableBonds', 'NumAromaticRings', 'TPSA', 'NumHeteroatoms',
        'NumRings', 'NumSaturatedRings', 'NumAliphaticRings', 'FractionCSP3',
        'MolarRefractivity', 'BertzCT', 'HeavyAtomCount'
    ]
    
    # Filter to only available columns
    feature_cols = [col for col in feature_cols if col in df.columns]
    print(f"Using {len(feature_cols)} features: {', '.join(feature_cols)}")
    
    # Initialize recommender
    recommender = CompoundRecommender(df, feature_cols, n_neighbors=10)
    
    # Generate batch recommendations
    recommendations = recommender.batch_recommend(n_alternatives=5)
    
    # Create report
    if not recommendations.empty:
        recommender.create_summary_report(recommendations)
    else:
        print("  WARNING: No recommendations generated")
    
    print("\n" + "="*80)
    print("[OK] RECOMMENDATION SYSTEM COMPLETE")
    print("="*80)
    print("\nOutputs:")
    print("  - outputs/analysis/alternatives.csv")
    print("  - outputs/ALTERNATIVES_REPORT.md")


