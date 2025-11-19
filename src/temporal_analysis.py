#!/usr/bin/env python3
"""
Temporal Analysis Module
========================

Analyze how pesticide toxicity to honey bees has changed over time (1832-2023).

Features:
- Temporal trend analysis using Mann-Kendall test
- Decade-by-decade comparison
- Rolling average calculations
- Pesticide type evolution over time
- Source comparison over time

Author: IME 372 Project Team
Date: November 2025
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import json
import os
from typing import Dict, Tuple
import warnings
warnings.filterwarnings('ignore')


class TemporalAnalyzer:
    """
    Analyze toxicity trends over time.
    """
    
    def __init__(self, df: pd.DataFrame):
        """
        Initialize temporal analyzer.
        
        Args:
            df: DataFrame with 'year' and 'label' columns
        """
        self.df = df.copy()
        
        # Validate required columns
        if 'year' not in self.df.columns or 'label' not in self.df.columns:
            raise ValueError("DataFrame must contain 'year' and 'label' columns")
    
    def mann_kendall_test(self, x: np.ndarray, y: np.ndarray) -> Tuple[float, float, str]:
        """
        Perform Mann-Kendall trend test.
        
        Args:
            x: Time values (years)
            y: Toxicity rates
            
        Returns:
            Tuple of (tau, p_value, trend_direction)
        """
        tau, p_value = stats.kendalltau(x, y)
        
        if p_value < 0.05:
            if tau > 0:
                trend_direction = "increasing (significant)"
            else:
                trend_direction = "decreasing (significant)"
        else:
            trend_direction = "no significant trend"
        
        return tau, p_value, trend_direction
    
    def toxicity_by_year(self, save_path='outputs/figures/', min_compounds=5) -> Dict:
        """
        Analyze how toxicity has changed by year.
        
        Args:
            save_path: Directory to save figures
            min_compounds: Minimum number of compounds per year to include
            
        Returns:
            Dictionary with analysis results
        """
        print("\n" + "="*80)
        print("TEMPORAL TREND ANALYSIS")
        print("="*80)
        
        # Create output directory
        os.makedirs(save_path, exist_ok=True)
        os.makedirs('outputs/analysis', exist_ok=True)
        
        # Calculate toxicity rate by year
        yearly = self.df.groupby('year').agg({
            'label': ['mean', 'count', 'sum']
        }).reset_index()
        yearly.columns = ['year', 'toxicity_rate', 'n_compounds', 'n_toxic']
        yearly['toxicity_pct'] = yearly['toxicity_rate'] * 100
        yearly['n_non_toxic'] = yearly['n_compounds'] - yearly['n_toxic']
        
        # Filter years with minimum compounds
        yearly = yearly[yearly['n_compounds'] >= min_compounds].copy()
        
        print(f"\nAnalyzing {len(yearly)} years with >={min_compounds} compounds")
        print(f"Year range: {yearly['year'].min()} - {yearly['year'].max()}")
        print(f"Total compounds: {yearly['n_compounds'].sum()}")
        
        # Calculate rolling average
        yearly['rolling_avg'] = yearly['toxicity_pct'].rolling(
            window=10, min_periods=3, center=True
        ).mean()
        
        # Mann-Kendall trend test
        tau, p_value, trend_direction = self.mann_kendall_test(
            yearly['year'].values, 
            yearly['toxicity_pct'].values
        )
        
        print(f"\n[DATA] Trend Analysis:")
        print(f"  Mann-Kendall  : {tau:.4f}")
        print(f"  p-value: {p_value:.4f}")
        print(f"  Trend: {trend_direction}")
        
        # Linear regression
        slope, intercept, r_value, p_val_reg, std_err = stats.linregress(
            yearly['year'], yearly['toxicity_pct']
        )
        
        print(f"\n[CHART] Linear Regression:")
        print(f"  Slope: {slope:.4f}% per year")
        print(f"  R : {r_value**2:.4f}")
        print(f"  p-value: {p_val_reg:.4f}")
        
        # Create visualization
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))
        
        # Top plot: Toxicity percentage over time
        ax1.scatter(yearly['year'], yearly['toxicity_pct'], 
                   s=yearly['n_compounds']*3, alpha=0.5, 
                   c=yearly['toxicity_pct'], cmap='RdYlGn_r',
                   label='Annual toxicity rate', edgecolors='black', linewidth=0.5)
        
        # Rolling average
        ax1.plot(yearly['year'], yearly['rolling_avg'], 
                color='red', linewidth=3, label='10-year rolling average',
                alpha=0.8)
        
        # Linear trend
        trend_line = intercept + slope * yearly['year']
        ax1.plot(yearly['year'], trend_line, 
                '--', color='black', linewidth=2, alpha=0.7,
                label=f'Linear trend (slope={slope:.3f}%/yr)')
        
        ax1.set_xlabel('Publication Year', fontsize=12, fontweight='bold')
        ax1.set_ylabel('% Toxic Compounds', fontsize=12, fontweight='bold')
        ax1.set_title(f'Temporal Trend in Pesticide Toxicity to Honey Bees (1832-2023)\n'
                     f'Trend: {trend_direction.capitalize()} (τ={tau:.3f}, p={p_value:.4f})',
                     fontsize=14, fontweight='bold')
        ax1.legend(loc='best', fontsize=10)
        ax1.grid(alpha=0.3, linestyle='--')
        ax1.set_ylim(0, 100)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap='RdYlGn_r', 
                                    norm=plt.Normalize(vmin=0, vmax=100))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax1, label='Toxicity %')
        
        # Bottom plot: Number of compounds over time
        ax2.bar(yearly['year'], yearly['n_non_toxic'], 
               label='Non-toxic', color='#2ecc71', alpha=0.8)
        ax2.bar(yearly['year'], yearly['n_toxic'], 
               bottom=yearly['n_non_toxic'],
               label='Toxic', color='#e74c3c', alpha=0.8)
        
        ax2.set_xlabel('Publication Year', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of Compounds', fontsize=12, fontweight='bold')
        ax2.set_title('Compound Count by Year', fontsize=13, fontweight='bold')
        ax2.legend(loc='best', fontsize=10)
        ax2.grid(alpha=0.3, linestyle='--', axis='y')
        
        plt.tight_layout()
        plt.savefig(f'{save_path}temporal_trend.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\nOK Figure saved: {save_path}temporal_trend.png")
        
        # Decade analysis
        self.df['decade'] = (self.df['year'] // 10) * 10
        decade_stats = self.df.groupby('decade').agg({
            'label': ['mean', 'count', 'sum']
        }).reset_index()
        decade_stats.columns = ['decade', 'toxicity_rate', 'n_compounds', 'n_toxic']
        decade_stats['toxicity_pct'] = decade_stats['toxicity_rate'] * 100
        
        print(f"\n  Decade Analysis:")
        print(decade_stats.to_string(index=False))
        
        # Save results
        results = {
            'mann_kendall_tau': float(tau),
            'mann_kendall_p_value': float(p_value),
            'trend_direction': trend_direction,
            'significant': bool(p_value < 0.05),
            'linear_slope': float(slope),
            'linear_r_squared': float(r_value**2),
            'years_analyzed': int(len(yearly)),
            'year_range': {
                'min': int(yearly['year'].min()),
                'max': int(yearly['year'].max())
            },
            'overall_toxicity_rate': float(self.df['label'].mean()),
            'yearly_data': yearly.to_dict('records'),
            'decade_data': decade_stats.to_dict('records')
        }
        
        # Save to JSON
        with open('outputs/analysis/temporal_trends.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\nOK Results saved: outputs/analysis/temporal_trends.json")
        
        return results
    
    def pesticide_type_evolution(self, save_path='outputs/figures/') -> Dict:
        """
        Analyze how different pesticide types evolved over time.
        
        Returns:
            Dictionary with evolution analysis for each pesticide type
        """
        print("\n" + "="*80)
        print("PESTICIDE TYPE EVOLUTION ANALYSIS")
        print("="*80)
        
        type_cols = ['herbicide', 'fungicide', 'insecticide']
        results = {}
        
        fig, axes = plt.subplots(len(type_cols), 1, figsize=(14, 12))
        
        for idx, ptype in enumerate(type_cols):
            # Filter to this type
            type_df = self.df[self.df[ptype] == 1].copy()
            
            if len(type_df) < 10:
                print(f"\nWARNING: Skipping {ptype}: insufficient data ({len(type_df)} compounds)")
                continue
            
            # Toxicity by year
            yearly = type_df.groupby('year').agg({
                'label': ['mean', 'count']
            }).reset_index()
            yearly.columns = ['year', 'toxicity_rate', 'n_compounds']
            yearly['toxicity_pct'] = yearly['toxicity_rate'] * 100
            
            # Filter years with at least 3 compounds
            yearly = yearly[yearly['n_compounds'] >= 3]
            
            if len(yearly) >= 10:
                # Trend test
                tau, p_value, trend = self.mann_kendall_test(
                    yearly['year'].values,
                    yearly['toxicity_pct'].values
                )
                
                results[ptype] = {
                    'tau': float(tau),
                    'p_value': float(p_value),
                    'trend': trend,
                    'n_years': int(len(yearly)),
                    'n_compounds': int(len(type_df)),
                    'mean_toxicity_pct': float(type_df['label'].mean() * 100)
                }
                
                print(f"\n{ptype.upper()}:")
                print(f"  Compounds: {len(type_df)}")
                print(f"  Toxicity rate: {results[ptype]['mean_toxicity_pct']:.1f}%")
                print(f"  Trend: {trend} ( ={tau:.3f}, p={p_value:.4f})")
                
                # Plot
                ax = axes[idx]
                ax.scatter(yearly['year'], yearly['toxicity_pct'],
                          s=yearly['n_compounds']*5, alpha=0.6,
                          label=f'{ptype.capitalize()} (n={len(type_df)})')
                
                # Rolling average
                if len(yearly) >= 5:
                    rolling = yearly['toxicity_pct'].rolling(window=5, min_periods=2, center=True).mean()
                    ax.plot(yearly['year'], rolling, linewidth=2, alpha=0.8)
                
                ax.set_ylabel('% Toxic', fontsize=11, fontweight='bold')
                ax.set_title(f'{ptype.capitalize()}: {trend}', fontsize=12)
                ax.legend(loc='best')
                ax.grid(alpha=0.3)
                ax.set_ylim(0, 100)
        
        axes[-1].set_xlabel('Publication Year', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{save_path}pesticide_type_evolution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\nOK Figure saved: {save_path}pesticide_type_evolution.png")
        
        return results


def run_temporal_analysis(data_path='outputs/dataset_with_descriptors.csv'):
    """
    Run complete temporal analysis pipeline.
    
    Args:
        data_path: Path to dataset CSV
    """
    print("="*80)
    print("COMPREHENSIVE TEMPORAL ANALYSIS")
    print("="*80)
    
    # Load data
    df = pd.read_csv(data_path)
    print(f"\nOK Loaded {len(df)} compounds from {data_path}")
    
    # Initialize analyzer
    analyzer = TemporalAnalyzer(df)
    
    # Run analyses
    yearly_results = analyzer.toxicity_by_year()
    type_results = analyzer.pesticide_type_evolution()
    
    # Save combined results
    combined_results = {
        'overall_trends': yearly_results,
        'pesticide_type_trends': type_results,
        'analysis_date': pd.Timestamp.now().isoformat()
    }
    
    with open('outputs/analysis/temporal_analysis_complete.json', 'w') as f:
        json.dump(combined_results, f, indent=2)
    
    print("\n" + "="*80)
    print("OK TEMPORAL ANALYSIS COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print("  - outputs/figures/temporal_trend.png")
    print("  - outputs/figures/pesticide_type_evolution.png")
    print("  - outputs/analysis/temporal_trends.json")
    print("  - outputs/analysis/temporal_analysis_complete.json")
    
    return combined_results


if __name__ == "__main__":
    run_temporal_analysis()


