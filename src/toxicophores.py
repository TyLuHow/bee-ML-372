#!/usr/bin/env python3
"""
Toxicophore Identification System
==================================

Identifies structural alerts (toxicophores) associated with bee toxicity
using SMARTS pattern matching and statistical correlation analysis.

Author: IME 372 Project Team
Date: November 2025
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from collections import Counter
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import json
import os
import warnings
warnings.filterwarnings('ignore')


class ToxicophoreAnalyzer:
    """
    Identify and analyze structural alerts (toxicophores) associated with bee toxicity.
    
    Uses SMARTS patterns to detect common pesticide functional groups and
    correlates their presence with toxicity outcomes.
    """
    
    def __init__(self):
        """Initialize with comprehensive pesticide toxicophore patterns."""
        
        # Core pesticide toxicophores with SMARTS patterns
        self.toxicophore_patterns = {
            # Insecticide classes
            'organophosphate': '[P](=O)([O,S])[O,S]',
            'carbamate': 'NC(=O)O',
            'pyrethroid_ester': 'C(=O)OC',
            'neonicotinoid_core': 'c1ncc([C,N])cn1',
            'phenylpyrazole': 'c1ccc(cc1)n2ncc(c2)',
            
            # Functional groups
            'nitro_group': '[N+](=O)[O-]',
            'cyano_group': 'C#N',
            'aromatic_halogen': 'c[F,Cl,Br,I]',
            'chlorinated_aromatic': 'c(Cl)',
            'fluorinated': 'C(F)(F)F',
            
            # Ring systems
            'triazole': 'c1ncnn1',
            'imidazole': 'c1nccn1',
            'pyridine': 'c1ccncc1',
            
            # Reactive groups
            'sulfonyl': 'S(=O)(=O)',
            'phosphate': 'P(=O)(O)(O)',
            'aromatic_amine': 'cN',
            'phenol': 'c[OH]',
            
            # Specific pesticide markers
            'urea_derivative': 'NC(=O)N',
            'oxime': 'C=NO',
            'methylenedioxyphenyl': 'c1ccc2OCOc2c1',
        }
        
        # Human-readable names
        self.toxicophore_names = {
            'organophosphate': 'Organophosphate',
            'carbamate': 'Carbamate',
            'pyrethroid_ester': 'Pyrethroid Ester',
            'neonicotinoid_core': 'Neonicotinoid',
            'phenylpyrazole': 'Phenylpyrazole',
            'nitro_group': 'Nitro Group',
            'cyano_group': 'Cyano Group',
            'aromatic_halogen': 'Aromatic Halogen',
            'chlorinated_aromatic': 'Chlorinated Aromatic',
            'fluorinated': 'Trifluoromethyl',
            'triazole': 'Triazole Ring',
            'imidazole': 'Imidazole Ring',
            'pyridine': 'Pyridine Ring',
            'sulfonyl': 'Sulfonyl',
            'phosphate': 'Phosphate',
            'aromatic_amine': 'Aromatic Amine',
            'phenol': 'Phenol',
            'urea_derivative': 'Urea Derivative',
            'oxime': 'Oxime',
            'methylenedioxyphenyl': 'Methylenedioxyphenyl'
        }
    
    def find_toxicophores(self, smiles: str) -> Dict[str, bool]:
        """
        Identify which toxicophores are present in a molecule.
        
        Args:
            smiles: SMILES string of molecule
            
        Returns:
            Dictionary mapping toxicophore names to presence (True/False)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {name: False for name in self.toxicophore_patterns.keys()}
        
        toxicophores = {}
        for name, smarts in self.toxicophore_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                toxicophores[name] = mol.HasSubstructMatch(pattern)
            else:
                toxicophores[name] = False
        
        return toxicophores
    
    def analyze_dataset(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Find toxicophores in all compounds and correlate with toxicity.
        
        Args:
            df: DataFrame with 'SMILES' and 'label' columns
            
        Returns:
            DataFrame with toxicophore correlation statistics
        """
        print("[ANALYSIS] Analyzing toxicophores in dataset...")
        
        # Find toxicophores for each compound
        toxicophore_data = []
        failed_smiles = 0
        
        for idx, row in df.iterrows():
            tox_dict = self.find_toxicophores(row['SMILES'])
            if tox_dict is None or all(v is False for v in tox_dict.values()):
                failed_smiles += 1
            
            tox_dict['label'] = row['label']
            tox_dict['name'] = row.get('name', f'Compound_{idx}')
            tox_dict['CID'] = row.get('CID', '')
            toxicophore_data.append(tox_dict)
        
        tox_df = pd.DataFrame(toxicophore_data)
        
        print(f"  OK Analyzed {len(df)} compounds")
        print(f"  WARNING: Failed to parse: {failed_smiles}")
        
        # Calculate statistics for each toxicophore
        correlations = []
        
        for col in self.toxicophore_patterns.keys():
            if col in tox_df.columns:
                # Contingency table
                ct = pd.crosstab(tox_df[col], tox_df['label'])
                
                if True in ct.index and False in ct.index:
                    # Calculate metrics
                    n_with = ct.loc[True].sum()
                    n_without = ct.loc[False].sum()
                    toxic_with = ct.loc[True, 1] if 1 in ct.columns else 0
                    toxic_without = ct.loc[False, 1] if 1 in ct.columns else 0
                    
                    toxicity_rate_with = toxic_with / n_with if n_with > 0 else 0
                    toxicity_rate_without = toxic_without / n_without if n_without > 0 else 0
                    
                    # Enrichment ratio
                    enrichment = (toxicity_rate_with / toxicity_rate_without 
                                 if toxicity_rate_without > 0 else float('inf'))
                    
                    # Chi-square test
                    chi2, p_value, _, _ = stats.chi2_contingency(ct)
                    
                    # Odds ratio
                    a, b = toxic_with, n_with - toxic_with
                    c, d = toxic_without, n_without - toxic_without
                    odds_ratio = (a * d) / (b * c) if (b * c) > 0 else float('inf')
                    
                    correlations.append({
                        'toxicophore': col,
                        'display_name': self.toxicophore_names[col],
                        'n_compounds': int(n_with),
                        'n_toxic': int(toxic_with),
                        'n_nontoxic': int(n_with - toxic_with),
                        'toxicity_rate': round(toxicity_rate_with * 100, 2),
                        'baseline_rate': round(toxicity_rate_without * 100, 2),
                        'enrichment': round(enrichment, 2) if enrichment != float('inf') else 999.0,
                        'odds_ratio': round(odds_ratio, 2) if odds_ratio != float('inf') else 999.0,
                        'chi_square': round(chi2, 4),
                        'p_value': round(p_value, 4),
                        'significant': p_value < 0.05
                    })
        
        # Create results DataFrame
        corr_df = pd.DataFrame(correlations)
        corr_df = corr_df.sort_values('enrichment', ascending=False)
        
        # Save results
        os.makedirs('outputs/analysis', exist_ok=True)
        corr_df.to_csv('outputs/analysis/toxicophore_analysis.csv', index=False)
        
        # Save detailed results
        results = {
            'summary': {
                'n_compounds': len(df),
                'n_toxicophores': len(self.toxicophore_patterns),
                'n_significant': int(corr_df['significant'].sum()),
                'failed_parsing': int(failed_smiles)
            },
            'toxicophores': corr_df.to_dict('records')
        }
        
        with open('outputs/analysis/toxicophore_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n  OK Results saved to outputs/analysis/")
        print(f"  [DATA] Found {len(corr_df)} toxicophores")
        print(f"  [STAT] {corr_df['significant'].sum()} statistically significant (p<0.05)")
        
        return corr_df
    
    def create_summary_plots(self, corr_df: pd.DataFrame):
        """
        Create comprehensive visualization of toxicophore analysis.
        
        Args:
            corr_df: DataFrame from analyze_dataset()
        """
        print("\n[DATA] Creating toxicophore visualizations...")
        
        os.makedirs('outputs/figures', exist_ok=True)
        
        # Sort by enrichment
        plot_df = corr_df.nlargest(15, 'enrichment')
        
        # Figure 1: Enrichment plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        colors = ['#e74c3c' if sig else '#95a5a6' for sig in plot_df['significant']]
        bars = ax.barh(plot_df['display_name'], plot_df['enrichment'], color=colors, alpha=0.7)
        
        ax.axvline(x=1, color='black', linestyle='--', linewidth=1, label='No enrichment')
        ax.set_xlabel('Toxicity Enrichment Ratio', fontsize=12, fontweight='bold')
        ax.set_ylabel('Toxicophore', fontsize=12, fontweight='bold')
        ax.set_title('Toxicophore Enrichment in Toxic vs Non-toxic Compounds', 
                     fontsize=14, fontweight='bold', pad=20)
        ax.legend(['Baseline (1.0)', 'Significant (p<0.05)', 'Not significant'], 
                 loc='lower right')
        ax.grid(axis='x', alpha=0.3)
        
        # Add value labels
        for bar, val in zip(bars, plot_df['enrichment']):
            if val < 900:  # Don't show infinity values
                ax.text(val + 0.1, bar.get_y() + bar.get_height()/2, 
                       f'{val:.2f}x', va='center', fontsize=9)
        
        plt.tight_layout()
        plt.savefig('outputs/figures/toxicophore_enrichment.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Figure 2: Prevalence vs Toxicity Rate
        fig, ax = plt.subplots(figsize=(12, 8))
        
        scatter = ax.scatter(
            plot_df['n_compounds'], 
            plot_df['toxicity_rate'],
            s=plot_df['enrichment'].clip(upper=10) * 100,
            c=plot_df['enrichment'].clip(upper=10),
            cmap='RdYlGn_r',
            alpha=0.6,
            edgecolors='black',
            linewidth=1
        )
        
        # Add labels for significant toxicophores
        for idx, row in plot_df[plot_df['significant']].iterrows():
            ax.annotate(
                row['display_name'],
                xy=(row['n_compounds'], row['toxicity_rate']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3)
            )
        
        ax.set_xlabel('Number of Compounds with Toxicophore', fontsize=12, fontweight='bold')
        ax.set_ylabel('Toxicity Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title('Toxicophore Prevalence vs Toxicity Rate\n(Size = Enrichment)', 
                     fontsize=14, fontweight='bold', pad=20)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Enrichment Ratio', fontsize=10)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('outputs/figures/toxicophore_prevalence.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("  OK Saved 2 toxicophore visualization plots")
    
    def generate_report(self, corr_df: pd.DataFrame, output_path: str = 'outputs/TOXICOPHORE_REPORT.md'):
        """
        Generate comprehensive markdown report.
        
        Args:
            corr_df: Results from analyze_dataset()
            output_path: Path to save report
        """
        significant = corr_df[corr_df['significant']]
        top_5 = corr_df.nlargest(5, 'enrichment')
        
        report = f"""# Toxicophore Analysis Report
**Generated**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Executive Summary

- **Total Compounds Analyzed**: {len(corr_df)}
- **Toxicophores Identified**: {len(corr_df)}
- **Statistically Significant**: {len(significant)} (p < 0.05)
- **Highest Enrichment**: {top_5.iloc[0]['display_name']} ({top_5.iloc[0]['enrichment']:.2f}x)

## Top 5 Most Enriched Toxicophores

| Rank | Toxicophore | Compounds | Toxicity Rate | Enrichment | p-value |
|------|-------------|-----------|---------------|------------|---------|
"""
        
        for i, (idx, row) in enumerate(top_5.iterrows(), 1):
            report += f"| {i} | {row['display_name']} | {row['n_compounds']} | {row['toxicity_rate']:.1f}% | {row['enrichment']:.2f}x | {row['p_value']:.4f} |\n"
        
        report += f"""

## Statistically Significant Toxicophores (p < 0.05)

{len(significant)} toxicophores showed statistically significant association with bee toxicity:

| Toxicophore | Compounds | Toxic | Toxicity Rate | Baseline | Enrichment | p-value |
|-------------|-----------|-------|---------------|----------|------------|---------|
"""
        
        for idx, row in significant.iterrows():
            enrich_str = f"{row['enrichment']:.2f}x" if row['enrichment'] < 900 else "High"
            report += f"| {row['display_name']} | {row['n_compounds']} | {row['n_toxic']} | {row['toxicity_rate']:.1f}% | {row['baseline_rate']:.1f}% | {enrich_str} | {row['p_value']:.4f} |\n"
        
        report += """

## Key Findings

### High-Risk Structural Features

"""
        
        high_risk = corr_df[(corr_df['enrichment'] > 1.5) & (corr_df['significant'])]
        for idx, row in high_risk.iterrows():
            enrich_str = f"{row['enrichment']:.1f}x" if row['enrichment'] < 900 else "Very high"
            report += f"- **{row['display_name']}**: {enrich_str} enrichment ({row['toxicity_rate']:.1f}% toxicity rate)\n"
        
        report += """

## Visualizations

Generated plots (see `outputs/figures/`):
1. `toxicophore_enrichment.png` - Enrichment ratios for top toxicophores
2. `toxicophore_prevalence.png` - Prevalence vs toxicity rate scatter

## Methodology

**SMARTS Pattern Matching**: Toxicophores identified using RDKit substructure matching with predefined SMARTS patterns representing common pesticide functional groups.

**Statistical Testing**: Chi-square tests for independence between toxicophore presence and toxicity. Enrichment calculated as ratio of toxicity rates.

**Significance Threshold**: p < 0.05 (uncorrected for multiple testing)

---

*Report generated by ApisTox Toxicophore Analyzer*
"""
        
        # Save report
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(report)
        
        print(f"  OK Report saved to {output_path}")


# CLI interface for standalone execution
if __name__ == "__main__":
    print("="*80)
    print("TOXICOPHORE ANALYSIS")
    print("="*80)
    
    # Load dataset
    df = pd.read_csv('outputs/dataset_with_descriptors.csv')
    print(f"\n[FILE] Loaded dataset: {len(df)} compounds")
    
    # Run analysis
    analyzer = ToxicophoreAnalyzer()
    results = analyzer.analyze_dataset(df)
    
    # Create visualizations
    analyzer.create_summary_plots(results)
    
    # Generate report
    analyzer.generate_report(results)
    
    print("\n" + "="*80)
    print("[OK] TOXICOPHORE ANALYSIS COMPLETE")
    print("="*80)
    print("\nOutputs:")
    print("  - outputs/analysis/toxicophore_analysis.csv")
    print("  - outputs/analysis/toxicophore_results.json")
    print("  - outputs/figures/toxicophore_enrichment.png")
    print("  - outputs/figures/toxicophore_prevalence.png")
    print("  - outputs/TOXICOPHORE_REPORT.md")


