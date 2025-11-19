#!/usr/bin/env python3
"""
Master Analysis Script
======================

Runs all analysis modules in sequence:
1. Source Comparison (ECOTOX vs PPDB)
2. Toxicophore Identification
3. Alternative Recommendations
4. Comprehensive Analysis (Temporal + Chemical Space)
5. Scaffold Split Testing

Author: IME 372 Project Team
Date: November 2025
"""

import os
import sys
import subprocess
from pathlib import Path
import time


def run_script(script_path, description):
    """Run a Python script and report results."""
    print("\n" + "="*80)
    print(f">>  {description}")
    print("="*80)
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=False,
            text=True,
            check=True
        )
        
        elapsed = time.time() - start_time
        print(f"\n  {description} COMPLETE ({elapsed:.1f}s)")
        return True
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"\n  {description} FAILED ({elapsed:.1f}s)")
        print(f"  Error: {e}")
        return False
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n  {description} ERROR ({elapsed:.1f}s)")
        print(f"  Error: {e}")
        return False


def main():
    """Run all analyses in sequence."""
    print("="*80)
    print("RUNNING ALL ANALYSES")
    print("="*80)
    print("\nThis will execute all analysis modules:")
    print("  1. Source Comparison (ECOTOX vs PPDB)")
    print("  2. Toxicophore Identification")
    print("  3. Alternative Compound Recommendations")
    print("  4. Comprehensive Analysis (Temporal + Chemical Space)")
    print("  5. Scaffold Split Testing (Optional)")
    print("\nEstimated time: 10-15 minutes")
    print("="*80)
    
    # Track results
    results = {}
    total_start = time.time()
    
    # 1. Source Comparison
    if Path('src/source_comparison.py').exists():
        results['source_comparison'] = run_script(
            'src/source_comparison.py',
            'Source Comparison (ECOTOX vs PPDB)'
        )
    else:
        print("\n    src/source_comparison.py not found - skipping")
        results['source_comparison'] = None
    
    # 2. Toxicophore Analysis
    if Path('src/toxicophores.py').exists():
        results['toxicophores'] = run_script(
            'src/toxicophores.py',
            'Toxicophore Identification'
        )
    else:
        print("\n    src/toxicophores.py not found - skipping")
        results['toxicophores'] = None
    
    # 3. Alternative Recommendations
    if Path('src/recommendations.py').exists():
        results['recommendations'] = run_script(
            'src/recommendations.py',
            'Alternative Compound Recommendations'
        )
    else:
        print("\n    src/recommendations.py not found - skipping")
        results['recommendations'] = None
    
    # 4. Comprehensive Analysis
    if Path('run_comprehensive_analysis.py').exists():
        results['comprehensive'] = run_script(
            'run_comprehensive_analysis.py',
            'Comprehensive Analysis (Temporal + Chemical Space)'
        )
    else:
        print("\n    run_comprehensive_analysis.py not found - skipping")
        results['comprehensive'] = None
    
    # 5. Scaffold Split Testing (Optional - may fail if no SMILES)
    if Path('test_scaffold_split.py').exists():
        print("\n" + "="*80)
        print("    Scaffold Split Testing (Optional)")
        print("="*80)
        print("This may fail if SMILES column is not available in the dataset")
        
        results['scaffold_split'] = run_script(
            'test_scaffold_split.py',
            'Scaffold Split Testing'
        )
    else:
        print("\n    test_scaffold_split.py not found - skipping")
        results['scaffold_split'] = None
    
    # Summary
    total_elapsed = time.time() - total_start
    
    print("\n" + "="*80)
    print("ALL ANALYSES COMPLETE")
    print("="*80)
    
    print(f"\nTotal time: {total_elapsed/60:.1f} minutes")
    print("\n  Summary:")
    
    for name, status in results.items():
        if status is True:
            print(f"    {name}")
        elif status is False:
            print(f"    {name} (failed)")
        else:
            print(f"     {name} (skipped)")
    
    # Check outputs
    print("\n  Generated Outputs:")
    
    output_files = [
        'outputs/figures/source_comparison.png',
        'outputs/figures/toxicophore_enrichment.png',
        'outputs/figures/toxicophore_prevalence.png',
        'outputs/figures/temporal_trends.png',
        'outputs/figures/chemical_space_pca.png',
        'outputs/figures/chemical_space_tsne.png',
        'outputs/analysis/source_comparison.json',
        'outputs/analysis/toxicophore_results.json',
        'outputs/analysis/alternatives.csv',
        'outputs/analysis/temporal_analysis.json',
        'outputs/TOXICOPHORE_REPORT.md',
        'outputs/ALTERNATIVES_REPORT.md',
        'outputs/TEMPORAL_ANALYSIS_REPORT.md'
    ]
    
    found_count = 0
    for fpath in output_files:
        if Path(fpath).exists():
            size = Path(fpath).stat().st_size / 1024  # KB
            print(f"    {fpath} ({size:.1f} KB)")
            found_count += 1
        else:
            print(f"    {fpath} (not found)")
    
    print(f"\nGenerated {found_count}/{len(output_files)} expected outputs")
    
    # Success criteria
    success_count = sum(1 for v in results.values() if v is True)
    total_run = sum(1 for v in results.values() if v is not None)
    
    if success_count == total_run:
        print("\n  ALL ANALYSES SUCCESSFUL!")
        return 0
    elif success_count > 0:
        print(f"\n    {success_count}/{total_run} analyses successful")
        return 1
    else:
        print("\n  NO ANALYSES COMPLETED SUCCESSFULLY")
        return 2


if __name__ == "__main__":
    sys.exit(main())


