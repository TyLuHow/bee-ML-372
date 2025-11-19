#!/usr/bin/env python3
"""
Chemical Space Visualization Module
====================================

Visualize the chemical space of pesticides using dimensionality reduction techniques.

Features:
- PCA (2D and 3D projections)
- t-SNE visualization
- Interactive Plotly plots
- Clustering analysis
- Colored by toxicity, pesticide type, year

Author: IME 372 Project Team
Date: November 2025
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import json
import os
from typing import Dict, List
import warnings
warnings.filterwarnings('ignore')

# Try to import plotly for interactive plots
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("WARNING: Plotly not available. Only static plots will be generated.")


class ChemicalSpaceVisualizer:
    """
    Visualize chemical space using dimensionality reduction.
    """
    
    def __init__(self, df: pd.DataFrame, feature_cols: List[str]):
        """
        Initialize visualizer.
        
        Args:
            df: DataFrame with molecular descriptors
            feature_cols: List of column names for molecular features
        """
        self.df = df.copy()
        self.feature_cols = feature_cols
        
        # Validate features exist
        missing_cols = set(feature_cols) - set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing columns: {missing_cols}")
        
        # Extract and scale features
        self.X = df[feature_cols].values
        self.scaler = StandardScaler()
        self.X_scaled = self.scaler.fit_transform(self.X)
        
        print(f"OK Initialized with {len(df)} compounds and {len(feature_cols)} features")
    
    def perform_pca(self, n_components=3) -> Dict:
        """
        Perform PCA dimensionality reduction.
        
        Args:
            n_components: Number of principal components
            
        Returns:
            Dictionary with PCA results
        """
        print(f"\nPerforming PCA with {n_components} components...")
        
        pca = PCA(n_components=n_components, random_state=42)
        X_pca = pca.fit_transform(self.X_scaled)
        
        # Calculate explained variance
        explained_var = pca.explained_variance_ratio_
        cumulative_var = np.cumsum(explained_var)
        
        print(f"OK PCA complete")
        print(f"  Explained variance by component:")
        for i, var in enumerate(explained_var):
            print(f"    PC{i+1}: {var*100:.2f}% (cumulative: {cumulative_var[i]*100:.2f}%)")
        
        # Add PCA coordinates to dataframe
        for i in range(n_components):
            self.df[f'PC{i+1}'] = X_pca[:, i]
        
        return {
            'X_pca': X_pca,
            'explained_variance': explained_var.tolist(),
            'cumulative_variance': cumulative_var.tolist(),
            'components': pca.components_.tolist()
        }
    
    def perform_tsne(self, perplexity=30, n_iter=1000) -> np.ndarray:
        """
        Perform t-SNE dimensionality reduction.
        
        Args:
            perplexity: t-SNE perplexity parameter
            n_iter: Number of iterations
            
        Returns:
            2D t-SNE coordinates
        """
        print(f"\nPerforming t-SNE (perplexity={perplexity}, n_iter={n_iter})...")
        
        tsne = TSNE(n_components=2, perplexity=perplexity, max_iter=n_iter,
                random_state=42, verbose=0)
        X_tsne = tsne.fit_transform(self.X_scaled)
        
        # Add to dataframe
        self.df['tSNE1'] = X_tsne[:, 0]
        self.df['tSNE2'] = X_tsne[:, 1]
        
        print(f"OK t-SNE complete")
        
        return X_tsne
    
    def perform_clustering(self, n_clusters=5, method='kmeans') -> np.ndarray:
        """
        Perform clustering in chemical space.
        
        Args:
            n_clusters: Number of clusters
            method: Clustering method ('kmeans')
            
        Returns:
            Cluster labels
        """
        print(f"\nPerforming {method} clustering with {n_clusters} clusters...")
        
        if method == 'kmeans':
            clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            labels = clusterer.fit_predict(self.X_scaled)
            
            self.df['cluster'] = labels
            
            # Calculate cluster statistics
            cluster_stats = []
            for cluster_id in range(n_clusters):
                cluster_mask = labels == cluster_id
                toxicity_rate = self.df.loc[cluster_mask, 'label'].mean()
                cluster_stats.append({
                    'cluster': cluster_id,
                    'n_compounds': int(cluster_mask.sum()),
                    'toxicity_rate': float(toxicity_rate),
                    'toxicity_pct': float(toxicity_rate * 100)
                })
            
            cluster_df = pd.DataFrame(cluster_stats).sort_values('toxicity_rate', ascending=False)
            print(f"\nOK Clustering complete")
            print("\nCluster toxicity rates:")
            print(cluster_df.to_string(index=False))
            
            return labels
        else:
            raise ValueError(f"Unknown clustering method: {method}")
    
    def plot_pca_2d(self, save_path='outputs/figures/chemical_space_pca_2d.png'):
        """Create 2D PCA plot colored by toxicity."""
        if 'PC1' not in self.df.columns:
            self.perform_pca(n_components=3)
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot 1: Colored by toxicity
        scatter1 = axes[0].scatter(
            self.df['PC1'], self.df['PC2'],
            c=self.df['label'], cmap='RdYlGn_r',
            s=50, alpha=0.6, edgecolors='black', linewidth=0.5
        )
        axes[0].set_xlabel(f'PC1 ({self.df["PC1"].std():.2f} std)', fontsize=12)
        axes[0].set_ylabel(f'PC2 ({self.df["PC2"].std():.2f} std)', fontsize=12)
        axes[0].set_title('Chemical Space: PCA (colored by toxicity)', fontsize=13, fontweight='bold')
        axes[0].grid(alpha=0.3)
        plt.colorbar(scatter1, ax=axes[0], label='Toxic (1) vs Non-toxic (0)')
        
        # Plot 2: Colored by pesticide type
        pesticide_type = []
        for _, row in self.df.iterrows():
            if row['insecticide'] == 1:
                pesticide_type.append('Insecticide')
            elif row['herbicide'] == 1:
                pesticide_type.append('Herbicide')
            elif row['fungicide'] == 1:
                pesticide_type.append('Fungicide')
            else:
                pesticide_type.append('Other')
        
        self.df['pesticide_type'] = pesticide_type
        
        for ptype, color in [('Insecticide', '#e74c3c'), ('Herbicide', '#3498db'), 
                              ('Fungicide', '#2ecc71'), ('Other', '#95a5a6')]:
            mask = self.df['pesticide_type'] == ptype
            axes[1].scatter(
                self.df.loc[mask, 'PC1'], self.df.loc[mask, 'PC2'],
                c=color, label=ptype, s=50, alpha=0.6,
                edgecolors='black', linewidth=0.5
            )
        
        axes[1].set_xlabel(f'PC1', fontsize=12)
        axes[1].set_ylabel(f'PC2', fontsize=12)
        axes[1].set_title('Chemical Space: PCA (colored by pesticide type)', fontsize=13, fontweight='bold')
        axes[1].legend(loc='best')
        axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"OK 2D PCA plot saved: {save_path}")
    
    def plot_tsne_2d(self, save_path='outputs/figures/chemical_space_tsne.png'):
        """Create 2D t-SNE plot."""
        if 'tSNE1' not in self.df.columns:
            self.perform_tsne()
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot 1: Colored by toxicity
        scatter1 = axes[0].scatter(
            self.df['tSNE1'], self.df['tSNE2'],
            c=self.df['label'], cmap='RdYlGn_r',
            s=50, alpha=0.6, edgecolors='black', linewidth=0.5
        )
        axes[0].set_xlabel('t-SNE Component 1', fontsize=12)
        axes[0].set_ylabel('t-SNE Component 2', fontsize=12)
        axes[0].set_title('Chemical Space: t-SNE (colored by toxicity)', fontsize=13, fontweight='bold')
        axes[0].grid(alpha=0.3)
        plt.colorbar(scatter1, ax=axes[0], label='Toxic (1) vs Non-toxic (0)')
        
        # Plot 2: Colored by pesticide type
        for ptype, color in [('Insecticide', '#e74c3c'), ('Herbicide', '#3498db'), 
                              ('Fungicide', '#2ecc71'), ('Other', '#95a5a6')]:
            mask = self.df['pesticide_type'] == ptype
            axes[1].scatter(
                self.df.loc[mask, 'tSNE1'], self.df.loc[mask, 'tSNE2'],
                c=color, label=ptype, s=50, alpha=0.6,
                edgecolors='black', linewidth=0.5
            )
        
        axes[1].set_xlabel('t-SNE Component 1', fontsize=12)
        axes[1].set_ylabel('t-SNE Component 2', fontsize=12)
        axes[1].set_title('Chemical Space: t-SNE (colored by pesticide type)', fontsize=13, fontweight='bold')
        axes[1].legend(loc='best')
        axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"OK t-SNE plot saved: {save_path}")
    
    def create_interactive_plots(self, save_dir='outputs/figures/'):
        """Create interactive Plotly visualizations."""
        if not PLOTLY_AVAILABLE:
            print("WARNING: Plotly not available. Skipping interactive plots.")
            return
        
        if 'PC1' not in self.df.columns:
            self.perform_pca(n_components=3)
        
        # Interactive 2D PCA
        fig_pca = px.scatter(
            self.df, x='PC1', y='PC2', color='label',
            hover_data=['name', 'year', 'pesticide_type'] if 'name' in self.df.columns else None,
            title='Interactive Chemical Space: PCA Projection',
            labels={'label': 'Toxicity (1=Toxic, 0=Non-toxic)'},
            color_continuous_scale='RdYlGn_r',
            width=900, height=600
        )
        fig_pca.write_html(f'{save_dir}chemical_space_pca_interactive.html')
        print(f"OK Interactive PCA saved: {save_dir}chemical_space_pca_interactive.html")
        
        # Interactive 3D PCA
        fig_3d = px.scatter_3d(
            self.df, x='PC1', y='PC2', z='PC3', color='label',
            hover_data=['name', 'year', 'pesticide_type'] if 'name' in self.df.columns else None,
            title='Interactive 3D Chemical Space',
            labels={'label': 'Toxicity'},
            color_continuous_scale='RdYlGn_r',
            width=900, height=700
        )
        fig_3d.write_html(f'{save_dir}chemical_space_3d_interactive.html')
        print(f"OK Interactive 3D PCA saved: {save_dir}chemical_space_3d_interactive.html")


def run_chemical_space_analysis(data_path='outputs/dataset_with_descriptors.csv'):
    """
    Run complete chemical space analysis.
    
    Args:
        data_path: Path to dataset CSV
    """
    print("="*80)
    print("CHEMICAL SPACE VISUALIZATION")
    print("="*80)
    
    # Load data
    df = pd.read_csv(data_path)
    print(f"\nOK Loaded {len(df)} compounds")
    
    # Define molecular feature columns
    feature_cols = [
        'MolecularWeight', 'LogP', 'NumHDonors', 'NumHAcceptors',
        'NumRotatableBonds', 'NumAromaticRings', 'TPSA', 'NumHeteroatoms',
        'NumRings', 'NumSaturatedRings', 'NumAliphaticRings', 'FractionCSP3',
        'MolarRefractivity', 'BertzCT', 'HeavyAtomCount'
    ]
    
    # Filter to available columns
    available_features = [col for col in feature_cols if col in df.columns]
    print(f"Using {len(available_features)} molecular features")
    
    # Initialize visualizer
    visualizer = ChemicalSpaceVisualizer(df, available_features)
    
    # Run PCA
    pca_results = visualizer.perform_pca(n_components=3)
    
    # Run t-SNE
    visualizer.perform_tsne(perplexity=30, n_iter=1000)
    
    # Perform clustering
    visualizer.perform_clustering(n_clusters=5)
    
    # Create static plots
    visualizer.plot_pca_2d()
    visualizer.plot_tsne_2d()
    
    # Create interactive plots if available
    visualizer.create_interactive_plots()
    
    # Save results
    results = {
        'pca_explained_variance': pca_results['explained_variance'],
        'pca_cumulative_variance': pca_results['cumulative_variance'],
        'n_compounds': len(df),
        'n_features': len(available_features),
        'features_used': available_features
    }
    
    with open('outputs/analysis/chemical_space_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "="*80)
    print("OK CHEMICAL SPACE ANALYSIS COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print("  - outputs/figures/chemical_space_pca_2d.png")
    print("  - outputs/figures/chemical_space_tsne.png")
    if PLOTLY_AVAILABLE:
        print("  - outputs/figures/chemical_space_pca_interactive.html")
        print("  - outputs/figures/chemical_space_3d_interactive.html")
    print("  - outputs/analysis/chemical_space_results.json")
    
    return results


if __name__ == "__main__":
    run_chemical_space_analysis()


