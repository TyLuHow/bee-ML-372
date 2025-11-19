#!/usr/bin/env python3
"""
Data Explorer Backend Module
=============================

FastAPI router providing 8 endpoints for comprehensive dataset exploration
BEFORE users make predictions.

Endpoints:
1. GET /api/explorer/overview - Dataset statistics
2. GET /api/explorer/molecular-diversity - Molecular descriptor distributions
3. GET /api/explorer/toxicity-by-class - Toxicity by chemical type
4. GET /api/explorer/temporal-trends - Toxicity trends over time
5. GET /api/explorer/chemical-space - PCA and t-SNE coordinates
6. GET /api/explorer/toxicophores - Toxicophore enrichment
7. GET /api/explorer/correlations - Feature correlations
8. GET /api/explorer/property-distributions - 2D property relationships

Author: IME 372 Project Team
Date: November 2025
"""

from fastapi import APIRouter, HTTPException
from typing import Dict, Any
import pandas as pd
import numpy as np
from datetime import datetime
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

# Import existing modules
from src.toxicophores import ToxicophoreAnalyzer

# Import caching and utility functions
from app.backend.cache import cached
from app.backend.utils import load_data, calculate_confidence_interval

# Create router
router = APIRouter()

# Dataset path
DATASET_PATH = "data/raw/dataset_with_descriptors.csv"
DATA_VERSION = "1.0.0"

# Molecular descriptor columns
DESCRIPTOR_COLS = [
    'MolecularWeight', 'LogP', 'NumHDonors', 'NumHAcceptors',
    'NumRotatableBonds', 'NumAromaticRings', 'TPSA', 'NumHeteroatoms',
    'NumRings', 'NumSaturatedRings', 'NumAliphaticRings', 'FractionCSP3',
    'MolarRefractivity', 'BertzCT', 'HeavyAtomCount'
]

# Global variable for dataset (loaded once)
_dataset_cache = None


def load_dataset() -> pd.DataFrame:
    """
    Load dataset with caching.

    Returns:
        DataFrame with dataset
    """
    global _dataset_cache

    if _dataset_cache is None:
        _dataset_cache = load_data(DATASET_PATH)

    return _dataset_cache


# ==================== ENDPOINT 1: Overview ====================

@router.get("/overview")
@cached(cache_key_prefix="overview", ttl_seconds=3600)
def get_overview() -> Dict[str, Any]:
    """
    Get dataset overview statistics.

    Returns comprehensive statistics about the dataset including:
    - Total compounds
    - Temporal range
    - Class distribution (toxic vs non-toxic)
    - Source distribution
    - Chemical type distribution
    - Exposure type distribution
    """
    try:
        df = load_dataset()

        # Determine chemical type for each compound
        def get_chemical_type(row):
            if row['insecticide'] == 1:
                return 'insecticide'
            elif row['herbicide'] == 1:
                return 'herbicide'
            elif row['fungicide'] == 1:
                return 'fungicide'
            else:
                return 'other'

        df['chemical_type'] = df.apply(get_chemical_type, axis=1)

        # Calculate statistics
        total_compounds = len(df)
        toxic_count = int(df['label'].sum())
        non_toxic_count = total_compounds - toxic_count

        # Temporal range
        min_year = int(df['year'].min())
        max_year = int(df['year'].max())
        span_years = max_year - min_year

        # Source distribution
        source_dist = df['source'].value_counts().to_dict()
        source_dist = {k: int(v) for k, v in source_dist.items()}

        # Chemical type distribution
        chem_type_dist = df['chemical_type'].value_counts().to_dict()
        chem_type_dist = {k: int(v) for k, v in chem_type_dist.items()}

        # Exposure type distribution
        exposure_dist = df['toxicity_type'].value_counts().to_dict()
        exposure_dist = {k: int(v) for k, v in exposure_dist.items()}

        return {
            "total_compounds": total_compounds,
            "temporal_range": {
                "min": min_year,
                "max": max_year,
                "span_years": span_years
            },
            "class_distribution": {
                "toxic": toxic_count,
                "non_toxic": non_toxic_count,
                "ratio": round(toxic_count / non_toxic_count, 3) if non_toxic_count > 0 else 0
            },
            "source_distribution": source_dist,
            "chemical_types": chem_type_dist,
            "exposure_types": exposure_dist,
            "timestamp": datetime.now().isoformat(),
            "data_version": DATA_VERSION
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing overview: {str(e)}")


# ==================== ENDPOINT 2: Molecular Diversity ====================

@router.get("/molecular-diversity")
@cached(cache_key_prefix="molecular_diversity", ttl_seconds=3600)
def get_molecular_diversity() -> Dict[str, Any]:
    """
    Get molecular descriptor distributions with toxicity overlay.

    Returns distributions (histograms) of key molecular descriptors
    separated by toxic vs non-toxic compounds.
    """
    try:
        df = load_dataset()

        descriptors_to_analyze = ['MolecularWeight', 'LogP', 'TPSA', 'NumRings', 'NumHDonors', 'NumHAcceptors']
        results = []

        for descriptor in descriptors_to_analyze:
            if descriptor not in df.columns:
                continue

            # Separate by toxicity
            toxic_data = df[df['label'] == 1][descriptor].dropna()
            non_toxic_data = df[df['label'] == 0][descriptor].dropna()

            # Create histograms with consistent bins
            all_data = df[descriptor].dropna()
            min_val, max_val = all_data.min(), all_data.max()

            # Create 20 bins
            bins = np.linspace(min_val, max_val, 21)

            # Compute histograms
            toxic_hist, _ = np.histogram(toxic_data, bins=bins)
            non_toxic_hist, _ = np.histogram(non_toxic_data, bins=bins)

            # Create bin labels (ranges)
            bin_ranges = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(len(bins)-1)]

            # Overall statistics
            stats_dict = {
                "mean": float(all_data.mean()),
                "std": float(all_data.std()),
                "min": float(all_data.min()),
                "max": float(all_data.max()),
                "median": float(all_data.median())
            }

            results.append({
                "name": descriptor,
                "toxic": {
                    "bins": bin_ranges,
                    "counts": toxic_hist.tolist()
                },
                "non_toxic": {
                    "bins": bin_ranges,
                    "counts": non_toxic_hist.tolist()
                },
                "stats": stats_dict
            })

        return {
            "descriptors": results,
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing molecular diversity: {str(e)}")


# ==================== ENDPOINT 3: Toxicity by Class ====================

@router.get("/toxicity-by-class")
@cached(cache_key_prefix="toxicity_by_class", ttl_seconds=3600)
def get_toxicity_by_class() -> Dict[str, Any]:
    """
    Get toxicity rates by chemical type with statistical tests.

    Returns toxicity rates for insecticides, herbicides, fungicides, and other
    with confidence intervals and chi-square tests.
    """
    try:
        df = load_dataset()

        # Determine chemical type
        def get_chemical_type(row):
            if row['insecticide'] == 1:
                return 'insecticide'
            elif row['herbicide'] == 1:
                return 'herbicide'
            elif row['fungicide'] == 1:
                return 'fungicide'
            else:
                return 'other'

        df['chemical_type'] = df.apply(get_chemical_type, axis=1)

        # Calculate statistics for each type
        results = []
        overall_toxicity_rate = df['label'].mean()

        for chem_type in ['insecticide', 'herbicide', 'fungicide', 'other']:
            subset = df[df['chemical_type'] == chem_type]
            total = len(subset)

            if total == 0:
                continue

            toxic = int(subset['label'].sum())
            toxicity_rate = toxic / total

            # Calculate 95% confidence interval using utility function
            ci_lower, ci_upper = calculate_confidence_interval(toxicity_rate, total, confidence=0.95)

            # Chi-square test against overall distribution
            # H0: toxicity rate same as overall
            observed = [toxic, total - toxic]
            expected_toxic = total * overall_toxicity_rate
            expected_non_toxic = total * (1 - overall_toxicity_rate)
            expected = [expected_toxic, expected_non_toxic]

            if expected_toxic > 0 and expected_non_toxic > 0:
                chi_square_stat = sum((o - e)**2 / e for o, e in zip(observed, expected))
                p_value = 1 - stats.chi2.cdf(chi_square_stat, df=1)
            else:
                p_value = 1.0

            results.append({
                "name": chem_type,
                "total": total,
                "toxic": toxic,
                "toxicity_rate": round(toxicity_rate, 4),
                "confidence_interval": [round(ci_lower, 4), round(ci_upper, 4)],
                "chi_square_p_value": round(p_value, 6)
            })

        return {
            "chemical_classes": results,
            "overall_toxicity_rate": round(overall_toxicity_rate, 4),
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing toxicity by class: {str(e)}")


# ==================== ENDPOINT 4: Temporal Trends ====================

@router.get("/temporal-trends")
@cached(cache_key_prefix="temporal_trends", ttl_seconds=3600)
def get_temporal_trends() -> Dict[str, Any]:
    """
    Get toxicity trends over time.

    Returns decade-level analysis, rolling averages, and Mann-Kendall trend test.
    """
    try:
        df = load_dataset()

        # Group by decade
        df['decade'] = (df['year'] // 10) * 10
        df['decade_label'] = df['decade'].apply(lambda x: f"{x}s")

        # Calculate decade statistics
        decade_stats = df.groupby('decade').agg({
            'label': ['count', 'sum', 'mean']
        }).reset_index()
        decade_stats.columns = ['decade', 'count', 'toxic_count', 'toxicity_rate']

        decades_result = []
        for _, row in decade_stats.iterrows():
            decade = int(row['decade'])
            decades_result.append({
                "decade": f"{decade}s",
                "year_start": decade,
                "year_end": decade + 9,
                "count": int(row['count']),
                "toxic_count": int(row['toxic_count']),
                "toxicity_rate": round(row['toxicity_rate'], 4)
            })

        # Calculate rolling average by year (for years with >= 5 compounds)
        yearly = df.groupby('year').agg({
            'label': ['mean', 'count']
        }).reset_index()
        yearly.columns = ['year', 'toxicity_rate', 'count']
        yearly = yearly[yearly['count'] >= 5].sort_values('year')

        if len(yearly) >= 10:
            # 10-year rolling average
            yearly['rolling_avg'] = yearly['toxicity_rate'].rolling(
                window=10, min_periods=3, center=True
            ).mean()

            rolling_years = yearly['year'].tolist()
            rolling_rates = yearly['rolling_avg'].fillna(yearly['toxicity_rate']).tolist()
        else:
            rolling_years = yearly['year'].tolist()
            rolling_rates = yearly['toxicity_rate'].tolist()

        # Mann-Kendall trend test
        if len(yearly) >= 10:
            tau, p_value = stats.kendalltau(yearly['year'], yearly['toxicity_rate'])

            if p_value < 0.05:
                if tau > 0:
                    trend = "increasing"
                else:
                    trend = "decreasing"
            else:
                trend = "no trend"
        else:
            tau = 0.0
            p_value = 1.0
            trend = "insufficient data"

        return {
            "decades": decades_result,
            "rolling_average": {
                "years": rolling_years,
                "rates": [round(r, 4) for r in rolling_rates]
            },
            "mann_kendall_test": {
                "tau": round(tau, 4),
                "p_value": round(p_value, 6),
                "trend": trend
            },
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing temporal trends: {str(e)}")


# ==================== ENDPOINT 5: Chemical Space ====================

@router.get("/chemical-space")
@cached(cache_key_prefix="chemical_space", ttl_seconds=3600)
def get_chemical_space() -> Dict[str, Any]:
    """
    Get PCA and t-SNE coordinates for chemical space visualization.

    Returns 2D PCA and t-SNE coordinates with compound metadata.
    """
    try:
        df = load_dataset()

        # Select descriptor columns
        feature_cols = [col for col in DESCRIPTOR_COLS if col in df.columns]

        # Extract features
        X = df[feature_cols].values

        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # PCA
        pca = PCA(n_components=2, random_state=42)
        X_pca = pca.fit_transform(X_scaled)

        # t-SNE (with lower perplexity for smaller datasets)
        perplexity = min(30, len(df) // 3)
        tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42, max_iter=1000)
        X_tsne = tsne.fit_transform(X_scaled)

        # Determine chemical type
        def get_chemical_type(row):
            if row['insecticide'] == 1:
                return 'insecticide'
            elif row['herbicide'] == 1:
                return 'herbicide'
            elif row['fungicide'] == 1:
                return 'fungicide'
            else:
                return 'other'

        df['chemical_type'] = df.apply(get_chemical_type, axis=1)

        # Build compound data
        compound_data = []
        for idx, row in df.iterrows():
            compound_data.append({
                "cid": int(row['CID']) if pd.notna(row['CID']) else 0,
                "name": str(row['name']) if 'name' in row and pd.notna(row['name']) else "",
                "toxic": bool(row['label']),
                "year": int(row['year']),
                "chemical_type": row['chemical_type'],
                "LogP": float(row['LogP']) if 'LogP' in row else 0.0
            })

        # PCA coordinates
        pca_coords = [[float(x), float(y)] for x, y in X_pca]

        # t-SNE coordinates
        tsne_coords = [[float(x), float(y)] for x, y in X_tsne]

        return {
            "pca": {
                "coordinates": pca_coords,
                "variance_explained": [round(v, 4) for v in pca.explained_variance_ratio_.tolist()],
                "compound_data": compound_data
            },
            "tsne": {
                "coordinates": tsne_coords,
                "perplexity": perplexity
            },
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing chemical space: {str(e)}")


# ==================== ENDPOINT 6: Toxicophores ====================

@router.get("/toxicophores")
@cached(cache_key_prefix="toxicophores", ttl_seconds=3600)
def get_toxicophores() -> Dict[str, Any]:
    """
    Get toxicophore enrichment data.

    Returns structural alerts enriched in toxic compounds with statistical tests.
    """
    try:
        df = load_dataset()

        # Initialize toxicophore analyzer
        analyzer = ToxicophoreAnalyzer()

        # Analyze toxicophores (this may take a while, hence caching)
        print("Analyzing toxicophores...")
        results = []

        for name, smarts in analyzer.toxicophore_patterns.items():
            from rdkit import Chem

            # Count presence in toxic vs non-toxic
            toxic_df = df[df['label'] == 1]
            non_toxic_df = df[df['label'] == 0]

            toxic_with = 0
            non_toxic_with = 0

            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue

            # Count in toxic compounds
            for smiles in toxic_df['SMILES']:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(pattern):
                    toxic_with += 1

            # Count in non-toxic compounds
            for smiles in non_toxic_df['SMILES']:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(pattern):
                    non_toxic_with += 1

            total_with = toxic_with + non_toxic_with
            total_without = len(df) - total_with

            if total_with == 0:
                continue

            # Calculate rates
            toxic_rate_with = toxic_with / total_with if total_with > 0 else 0
            toxic_rate_without = (len(toxic_df) - toxic_with) / total_without if total_without > 0 else 0

            # Enrichment ratio
            enrichment = toxic_rate_with / toxic_rate_without if toxic_rate_without > 0 else 999.0

            # Chi-square test
            if total_with > 5 and total_without > 5:
                contingency = [
                    [toxic_with, non_toxic_with],
                    [len(toxic_df) - toxic_with, len(non_toxic_df) - non_toxic_with]
                ]
                chi2, p_value, _, _ = stats.chi2_contingency(contingency)
            else:
                p_value = 1.0

            # Calculate confidence interval using utility function
            ci_lower, ci_upper = calculate_confidence_interval(toxic_rate_with, total_with, confidence=0.95)

            results.append({
                "name": analyzer.toxicophore_names[name],
                "smarts": smarts,
                "prevalence": round(total_with / len(df), 4),
                "toxic_rate_with": round(toxic_rate_with, 4),
                "toxic_rate_without": round(toxic_rate_without, 4),
                "enrichment_ratio": round(enrichment, 4) if enrichment < 999 else 999.0,
                "p_value": round(p_value, 6),
                "confidence_interval": [round(ci_lower, 4), round(ci_upper, 4)]
            })

        # Sort by enrichment and take top 10
        results.sort(key=lambda x: x['enrichment_ratio'], reverse=True)
        results = results[:10]

        return {
            "toxicophores": results,
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing toxicophores: {str(e)}")


# ==================== ENDPOINT 7: Correlations ====================

@router.get("/correlations")
@cached(cache_key_prefix="correlations", ttl_seconds=3600)
def get_correlations() -> Dict[str, Any]:
    """
    Get feature correlation matrix.

    Returns correlation matrix, top correlations, and network data for visualization.
    """
    try:
        df = load_dataset()

        # Select descriptor columns
        feature_cols = [col for col in DESCRIPTOR_COLS if col in df.columns]

        # Calculate correlation matrix
        corr_matrix = df[feature_cols].corr()

        # Convert to list format
        corr_values = corr_matrix.values.tolist()

        # Find top correlations (excluding diagonal)
        top_correlations = []
        for i in range(len(feature_cols)):
            for j in range(i+1, len(feature_cols)):
                corr_val = corr_matrix.iloc[i, j]
                top_correlations.append({
                    "feature1": feature_cols[i],
                    "feature2": feature_cols[j],
                    "correlation": round(corr_val, 4)
                })

        # Sort by absolute correlation
        top_correlations.sort(key=lambda x: abs(x['correlation']), reverse=True)
        top_correlations = top_correlations[:20]

        # Create network data (edges with |r| > 0.5)
        nodes = [{"id": feat, "label": feat} for feat in feature_cols]
        edges = []

        for i in range(len(feature_cols)):
            for j in range(i+1, len(feature_cols)):
                corr_val = corr_matrix.iloc[i, j]
                if abs(corr_val) > 0.5:
                    edges.append({
                        "source": feature_cols[i],
                        "target": feature_cols[j],
                        "weight": round(corr_val, 4)
                    })

        return {
            "correlation_matrix": {
                "features": feature_cols,
                "values": [[round(v, 4) for v in row] for row in corr_values]
            },
            "top_correlations": top_correlations,
            "network_data": {
                "nodes": nodes,
                "edges": edges
            },
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing correlations: {str(e)}")


# ==================== ENDPOINT 8: Property Distributions ====================

@router.get("/property-distributions")
@cached(cache_key_prefix="property_distributions", ttl_seconds=3600)
def get_property_distributions() -> Dict[str, Any]:
    """
    Get 2D property relationship scatter plots.

    Returns scatter plot data for key property relationships.
    """
    try:
        df = load_dataset()

        # Define property pairs to analyze
        property_pairs = [
            ('LogP', 'MolecularWeight'),
            ('TPSA', 'NumHDonors'),
            ('NumRotatableBonds', 'FractionCSP3')
        ]

        scatter_data = []

        for x_prop, y_prop in property_pairs:
            if x_prop not in df.columns or y_prop not in df.columns:
                continue

            # Extract data
            points = []
            for idx, row in df.iterrows():
                if pd.notna(row[x_prop]) and pd.notna(row[y_prop]):
                    points.append({
                        "x": float(row[x_prop]),
                        "y": float(row[y_prop]),
                        "toxic": bool(row['label']),
                        "name": str(row['name']) if 'name' in row and pd.notna(row['name']) else "",
                        "cid": int(row['CID']) if pd.notna(row['CID']) else 0
                    })

            scatter_data.append({
                "name": f"{x_prop}_vs_{y_prop}",
                "x_label": x_prop,
                "y_label": y_prop,
                "points": points
            })

        return {
            "scatter_data": scatter_data,
            "timestamp": datetime.now().isoformat()
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing property distributions: {str(e)}")


# ==================== CACHE MANAGEMENT ENDPOINT ====================

@router.get("/cache/stats")
def get_cache_stats() -> Dict[str, Any]:
    """
    Get cache statistics (for debugging/monitoring).

    Returns:
        Cache size and other metrics
    """
    from app.backend.cache import get_cache_stats
    return get_cache_stats()


@router.post("/cache/clear")
def clear_cache() -> Dict[str, str]:
    """
    Clear all cached data (admin endpoint).

    Returns:
        Success message
    """
    from app.backend.cache import clear_cache
    clear_cache()
    return {"message": "Cache cleared successfully"}
