#!/usr/bin/env python3
"""
Comprehensive Analysis Runner
==============================

Runs all Phase 2 exploratory analyses:
1. Temporal trend analysis
2. Chemical space visualization
3. Source comparison

Author: IME 372 Project Team
Date: November 2025
"""

import sys
import os
from datetime import datetime

def run_all_analyses():
    """Run all exploratory analyses."""
    
    print("="*80)
    print("COMPREHENSIVE ANALYSIS SUITE")
    print("="*80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Import analysis modules
    try:
        from src.temporal_analysis import run_temporal_analysis
        from src.chemical_space import run_chemical_space_analysis
    except ImportError as e:
        print(f"X Import error: {e}")
        print("  Make sure all dependencies are installed.")
        return False
    
    # Run temporal analysis
    print("\n" + "="*80)
    print("PHASE 1: TEMPORAL TREND ANALYSIS")
    print("="*80)
    try:
        temporal_results = run_temporal_analysis()
        print("OK Temporal analysis complete")
    except Exception as e:
        print(f"X Temporal analysis failed: {e}")
        import traceback
        traceback.print_exc()
        temporal_results = None
    
    # Run chemical space visualization
    print("\n" + "="*80)
    print("PHASE 2: CHEMICAL SPACE VISUALIZATION")
    print("="*80)
    try:
        chemical_space_results = run_chemical_space_analysis()
        print("OK Chemical space analysis complete")
    except Exception as e:
        print(f"X Chemical space analysis failed: {e}")
        import traceback
        traceback.print_exc()
        chemical_space_results = None
    
    # Summary
    print("\n" + "="*80)
    print("ANALYSIS SUITE COMPLETE")
    print("="*80)
    
    success_count = sum([
        temporal_results is not None,
        chemical_space_results is not None
    ])
    
    print(f"\nOK {success_count}/2 analyses completed successfully")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print("\nGenerated outputs:")
    print("  Temporal Analysis:")
    print("    - outputs/figures/temporal_trend.png")
    print("    - outputs/figures/pesticide_type_evolution.png")
    print("    - outputs/analysis/temporal_trends.json")
    print("\n  Chemical Space:")
    print("    - outputs/figures/chemical_space_pca_2d.png")
    print("    - outputs/figures/chemical_space_tsne.png")
    print("    - outputs/figures/chemical_space_*_interactive.html (if Plotly available)")
    print("    - outputs/analysis/chemical_space_results.json")
    
    print("\nNext steps:")
    print("  1. Review generated visualizations in outputs/figures/")
    print("  2. Check analysis results in outputs/analysis/")
    print("  3. Run additional analyses as needed")
    
    return success_count == 2


if __name__ == "__main__":
    success = run_all_analyses()
    sys.exit(0 if success else 1)




