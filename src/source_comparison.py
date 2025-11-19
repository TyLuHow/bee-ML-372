#!/usr/bin/env python3
"""
ECOTOX vs PPDB Source Comparison Analysis
==========================================

Compare toxicity assessments between different data sources.

Features:
- Agreement rate for overlapping compounds
- Distribution differences (Chi-square test)
- Source-specific biases
- Temporal coverage comparison
- Property distribution comparison

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
from typing import Dict
import warnings
warnings.filterwarnings('ignore')


def compare_data_sources(df: pd.DataFrame) -> Dict:
    """
    Compare toxicity assessments between ECOTOX and PPDB sources.
    
    Analysis:
    1. Agreement rate for overlapping compounds (if any)
    2. Distribution differences (Chi-square test)
    3. Source-specific biases
    4. Temporal coverage comparison
    5. Property distribution comparison
    
    Args:
        df: DataFrame with 'source', 'label', 'year', molecular descriptors
        
    Returns:
        Dictionary with comparison results
    """
    print("="*80)
    print("ECOTOX vs PPDB SOURCE COMPARISON")
    print("="*80)
    
    # Overall source distribution
    print("\n[DATA] Source Distribution:")
    source_counts = df['source'].value_counts()
    print(source_counts)
    
    source_stats = df.groupby('source').agg({
        'label': ['mean', 'count', 'sum'],
        'year': ['min', 'max', 'mean', 'std']
    }).round(2)
    
    print("\n[DATA] Source Statistics:")
    print(source_stats)
    
    # Chi-square test for source vs toxicity
    contingency = pd.crosstab(df['source'], df['label'])
    print("\n  Contingency Table:")
    print(contingency)
    
    chi2, p_value, dof, expected = stats.chi2_contingency(contingency)
    
    print(f"\n[ANALYSIS] Chi-Square Test (Independence):")
    print(f"     statistic: {chi2:.4f}")
    print(f"  p-value: {p_value:.4f}")
    print(f"  Degrees of freedom: {dof}")
    print(f"  Significant difference: {'YES (p<0.05)' if p_value < 0.05 else 'NO (p>=0.05)'}")
    
    # Check for overlapping compounds (same CID in multiple sources)
    duplicate_cids = df[df.duplicated('CID', keep=False)]['CID'].unique()
    
    agreement_rate = None
    disagreement_examples = []
    
    if len(duplicate_cids) > 0:
        print(f"\n[PROCESS] Overlapping Compounds:")
        print(f"  Found {len(duplicate_cids)} compounds in multiple sources")
        
        comparison = []
        for cid in duplicate_cids:
            entries = df[df['CID'] == cid]
            if len(entries) >= 2:
                sources_list = entries['source'].tolist()
                labels_list = entries['label'].tolist()
                
                # Check all pairs
                for i in range(len(entries)):
                    for j in range(i+1, len(entries)):
                        comparison.append({
                            'CID': cid,
                            'name': entries.iloc[i]['name'],
                            'source_1': sources_list[i],
                            'source_2': sources_list[j],
                            'label_1': labels_list[i],
                            'label_2': labels_list[j],
                            'agreement': labels_list[i] == labels_list[j]
                        })
        
        if comparison:
            comparison_df = pd.DataFrame(comparison)
            agreement_rate = comparison_df['agreement'].mean() * 100
            
            print(f"  Agreement rate: {agreement_rate:.1f}%")
            print(f"  Agreements: {comparison_df['agreement'].sum()}")
            print(f"  Disagreements: {(~comparison_df['agreement']).sum()}")
            
            # Show disagreements
            disagreements = comparison_df[~comparison_df['agreement']]
            if len(disagreements) > 0:
                print(f"\n  Example disagreements:")
                for idx, row in disagreements.head(5).iterrows():
                    print(f"    - {row['name']} (CID: {row['CID']})")
                    print(f"      {row['source_1']}: {'Toxic' if row['label_1']==1 else 'Non-toxic'}")
                    print(f"      {row['source_2']}: {'Toxic' if row['label_2']==1 else 'Non-toxic'}")
                
                disagreement_examples = disagreements.head(10).to_dict('records')
    else:
        print(f"\n[PROCESS] No overlapping compounds found (all CIDs unique to one source)")
    
    # Create visualizations
    os.makedirs('outputs/figures', exist_ok=True)
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # 1. Toxicity distribution by source
    ax1 = fig.add_subplot(gs[0, 0])
    source_counts = pd.crosstab(df['source'], df['label'])
    source_pcts = source_counts.div(source_counts.sum(axis=1), axis=0) * 100
    
    x = np.arange(len(source_pcts))
    width = 0.35
    
    ax1.bar(x - width/2, source_pcts[0], width, label='Non-toxic', color='#2ecc71', alpha=0.8)
    ax1.bar(x + width/2, source_pcts[1], width, label='Toxic', color='#e74c3c', alpha=0.8)
    
    ax1.set_xlabel('Data Source', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Percentage (%)', fontsize=11, fontweight='bold')
    ax1.set_title('Toxicity Distribution by Source', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(source_pcts.index, rotation=0)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_ylim(0, 100)
    
    # 2. Temporal coverage
    ax2 = fig.add_subplot(gs[0, 1])
    for source in df['source'].unique():
        source_data = df[df['source'] == source]
        ax2.hist(source_data['year'], bins=40, alpha=0.6, label=source, edgecolor='black', linewidth=0.5)
    
    ax2.set_xlabel('Publication Year', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of Compounds', fontsize=11, fontweight='bold')
    ax2.set_title('Temporal Coverage by Source', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Pesticide type distribution
    ax3 = fig.add_subplot(gs[1, 0])
    type_cols = ['herbicide', 'fungicide', 'insecticide', 'other_agrochemical']
    type_data = []
    
    for source in df['source'].unique():
        source_df = df[df['source'] == source]
        for ptype in type_cols:
            if ptype in source_df.columns:
                count = source_df[ptype].sum()
                type_data.append({
                    'Source': source,
                    'Type': ptype.replace('_', ' ').capitalize(),
                    'Count': count
                })
    
    if type_data:
        type_df = pd.DataFrame(type_data)
        type_pivot = type_df.pivot(index='Source', columns='Type', values='Count')
        type_pivot.plot(kind='bar', ax=ax3, width=0.8)
        ax3.set_xlabel('Data Source', fontsize=11, fontweight='bold')
        ax3.set_ylabel('Number of Compounds', fontsize=11, fontweight='bold')
        ax3.set_title('Pesticide Type Distribution by Source', fontsize=12, fontweight='bold')
        ax3.set_xticklabels(ax3.get_xticklabels(), rotation=0)
        ax3.legend(title='Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax3.grid(axis='y', alpha=0.3)
    
    # 4. Molecular weight distribution
    ax4 = fig.add_subplot(gs[1, 1])
    for source in df['source'].unique():
        source_data = df[df['source'] == source]
        if 'MolecularWeight' in source_data.columns:
            ax4.hist(source_data['MolecularWeight'], bins=30, alpha=0.6, label=source, 
                    edgecolor='black', linewidth=0.5)
    
    ax4.set_xlabel('Molecular Weight (g/mol)', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax4.set_title('Molecular Weight Distribution', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(axis='y', alpha=0.3)
    
    # 5. Toxicity rate over time by source
    ax5 = fig.add_subplot(gs[2, :])
    for source in df['source'].unique():
        source_df = df[df['source'] == source]
        
        # Group by decade
        source_df['decade'] = (source_df['year'] // 10) * 10
        decade_tox = source_df.groupby('decade')['label'].agg(['mean', 'count'])
        decade_tox = decade_tox[decade_tox['count'] >= 3]  # At least 3 compounds
        
        if len(decade_tox) > 0:
            ax5.plot(decade_tox.index, decade_tox['mean'] * 100, 
                    marker='o', label=source, linewidth=2, markersize=6)
    
    ax5.set_xlabel('Decade', fontsize=11, fontweight='bold')
    ax5.set_ylabel('Toxicity Rate (%)', fontsize=11, fontweight='bold')
    ax5.set_title('Toxicity Rate Trends by Source Over Time', fontsize=12, fontweight='bold')
    ax5.legend()
    ax5.grid(alpha=0.3)
    ax5.set_ylim(0, 100)
    
    plt.suptitle('Data Source Comparison: ECOTOX vs PPDB', 
                fontsize=16, fontweight='bold', y=0.995)
    
    plt.savefig('outputs/figures/source_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nOK Visualization saved: outputs/figures/source_comparison.png")
    
    # Calculate additional statistics
    toxicity_by_source = df.groupby('source')['label'].agg(['mean', 'count', 'sum'])
    toxicity_by_source['toxicity_pct'] = toxicity_by_source['mean'] * 100
    
    # Save results
    results = {
        'summary': {
            'n_total': len(df),
            'sources': df['source'].unique().tolist(),
            'source_counts': source_counts.to_dict(),
            'chi_square_statistic': float(chi2),
            'chi_square_p_value': float(p_value),
            'significant_difference': bool(p_value < 0.05),
            'n_overlapping_compounds': int(len(duplicate_cids)),
            'agreement_rate_pct': float(agreement_rate) if agreement_rate is not None else None
        },
        'source_statistics': {
            source: {
                'n_compounds': int(row['count']),
                'n_toxic': int(row['sum']),
                'toxicity_rate_pct': float(row['toxicity_pct']),
                'year_range': {
                    'min': int(df[df['source']==source]['year'].min()),
                    'max': int(df[df['source']==source]['year'].max()),
                    'mean': float(df[df['source']==source]['year'].mean())
                }
            }
            for source, row in toxicity_by_source.iterrows()
        },
        'disagreement_examples': disagreement_examples if disagreement_examples else []
    }
    
    os.makedirs('outputs/analysis', exist_ok=True)
    with open('outputs/analysis/source_comparison.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"OK Results saved: outputs/analysis/source_comparison.json")
    
    return results


if __name__ == "__main__":
    print("="*80)
    print("SOURCE COMPARISON ANALYSIS")
    print("="*80)
    
    # Load dataset
    df = pd.read_csv('outputs/dataset_with_descriptors.csv')
    print(f"\n[FILE] Loaded {len(df)} compounds")
    
    # Run comparison
    results = compare_data_sources(df)
    
    print("\n" + "="*80)
    print("[OK] SOURCE COMPARISON COMPLETE")
    print("="*80)
    print("\nKey Findings:")
    print(f"  - Sources compared: {len(results['summary']['sources'])}")
    print(f"  - Statistical difference: {'YES' if results['summary']['significant_difference'] else 'NO'}")
    if results['summary']['agreement_rate_pct']:
        print(f"  - Agreement rate: {results['summary']['agreement_rate_pct']:.1f}%")
    print("\nOutputs:")
    print("  - outputs/figures/source_comparison.png")
    print("  - outputs/analysis/source_comparison.json")


