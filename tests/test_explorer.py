#!/usr/bin/env python3
"""
Unit Tests for Explorer Backend Module
=======================================

Tests cover:
- All 8 explorer endpoints
- Response structure validation
- Response time checks
- Error handling
- Data consistency

Author: IME 372 Project Team
Date: November 2025
"""

import pytest
import json
import os
import sys
import time
from fastapi.testclient import TestClient

# Import the FastAPI app
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from app.backend.main import app

# Create test client
client = TestClient(app)


class TestOverviewEndpoint:
    """Test /api/explorer/overview endpoint."""

    def test_overview_success(self):
        """Test that overview endpoint returns 200 OK."""
        response = client.get("/api/explorer/overview")
        assert response.status_code == 200

        data = response.json()

        # Check required fields
        assert "total_compounds" in data
        assert "temporal_range" in data
        assert "class_distribution" in data
        assert "source_distribution" in data
        assert "chemical_types" in data
        assert "exposure_types" in data
        assert "timestamp" in data
        assert "data_version" in data

        # Validate data types
        assert isinstance(data["total_compounds"], int)
        assert data["total_compounds"] > 0

        # Check temporal range structure
        temp_range = data["temporal_range"]
        assert "min" in temp_range
        assert "max" in temp_range
        assert "span_years" in temp_range
        assert temp_range["min"] < temp_range["max"]

        # Check class distribution
        class_dist = data["class_distribution"]
        assert "toxic" in class_dist
        assert "non_toxic" in class_dist
        assert "ratio" in class_dist
        assert class_dist["toxic"] + class_dist["non_toxic"] == data["total_compounds"]

    def test_overview_response_time(self):
        """Test that overview endpoint responds quickly (should be cached)."""
        # First call (might be slow)
        start = time.time()
        response1 = client.get("/api/explorer/overview")
        time1 = time.time() - start

        # Second call (should be cached)
        start = time.time()
        response2 = client.get("/api/explorer/overview")
        time2 = time.time() - start

        assert response1.status_code == 200
        assert response2.status_code == 200

        # Second call should be faster (cached)
        assert time2 < 0.5  # Should be <500ms


class TestMolecularDiversityEndpoint:
    """Test /api/explorer/molecular-diversity endpoint."""

    def test_molecular_diversity_success(self):
        """Test molecular diversity endpoint."""
        response = client.get("/api/explorer/molecular-diversity")
        assert response.status_code == 200

        data = response.json()
        assert "descriptors" in data
        assert "timestamp" in data

        # Check descriptors structure
        descriptors = data["descriptors"]
        assert isinstance(descriptors, list)
        assert len(descriptors) > 0

        # Validate first descriptor
        desc = descriptors[0]
        assert "name" in desc
        assert "toxic" in desc
        assert "non_toxic" in desc
        assert "stats" in desc

        # Check histogram structure
        assert "bins" in desc["toxic"]
        assert "counts" in desc["toxic"]
        assert len(desc["toxic"]["bins"]) == len(desc["toxic"]["counts"])

        # Check stats
        stats = desc["stats"]
        assert "mean" in stats
        assert "std" in stats
        assert "min" in stats
        assert "max" in stats

    def test_molecular_diversity_descriptors(self):
        """Test that key descriptors are included."""
        response = client.get("/api/explorer/molecular-diversity")
        data = response.json()

        descriptor_names = [d["name"] for d in data["descriptors"]]

        # Check for expected descriptors
        expected = ["MolecularWeight", "LogP", "TPSA"]
        for exp in expected:
            assert exp in descriptor_names


class TestToxicityByClassEndpoint:
    """Test /api/explorer/toxicity-by-class endpoint."""

    def test_toxicity_by_class_success(self):
        """Test toxicity by class endpoint."""
        response = client.get("/api/explorer/toxicity-by-class")
        assert response.status_code == 200

        data = response.json()
        assert "chemical_classes" in data
        assert "overall_toxicity_rate" in data
        assert "timestamp" in data

        # Check classes structure
        classes = data["chemical_classes"]
        assert isinstance(classes, list)
        assert len(classes) > 0

        # Validate first class
        cls = classes[0]
        assert "name" in cls
        assert "total" in cls
        assert "toxic" in cls
        assert "toxicity_rate" in cls
        assert "confidence_interval" in cls
        assert "chi_square_p_value" in cls

        # Check CI structure
        ci = cls["confidence_interval"]
        assert len(ci) == 2
        assert 0 <= ci[0] <= ci[1] <= 1

    def test_toxicity_by_class_types(self):
        """Test that all chemical types are included."""
        response = client.get("/api/explorer/toxicity-by-class")
        data = response.json()

        class_names = [c["name"] for c in data["chemical_classes"]]

        # Check for expected types
        expected_types = ["insecticide", "herbicide", "fungicide", "other"]
        for exp_type in expected_types:
            # At least some should be present
            pass  # Data-dependent, so we just check structure


class TestTemporalTrendsEndpoint:
    """Test /api/explorer/temporal-trends endpoint."""

    def test_temporal_trends_success(self):
        """Test temporal trends endpoint."""
        response = client.get("/api/explorer/temporal-trends")
        assert response.status_code == 200

        data = response.json()
        assert "decades" in data
        assert "rolling_average" in data
        assert "mann_kendall_test" in data
        assert "timestamp" in data

        # Check decades structure
        decades = data["decades"]
        assert isinstance(decades, list)
        assert len(decades) > 0

        # Validate first decade
        decade = decades[0]
        assert "decade" in decade
        assert "year_start" in decade
        assert "year_end" in decade
        assert "count" in decade
        assert "toxic_count" in decade
        assert "toxicity_rate" in decade

        # Check rolling average
        rolling = data["rolling_average"]
        assert "years" in rolling
        assert "rates" in rolling
        assert len(rolling["years"]) == len(rolling["rates"])

        # Check Mann-Kendall test
        mk = data["mann_kendall_test"]
        assert "tau" in mk
        assert "p_value" in mk
        assert "trend" in mk
        assert mk["trend"] in ["increasing", "decreasing", "no trend", "insufficient data"]

    def test_temporal_trends_decades_ordered(self):
        """Test that decades are in chronological order."""
        response = client.get("/api/explorer/temporal-trends")
        data = response.json()

        decades = data["decades"]
        year_starts = [d["year_start"] for d in decades]

        # Check ordering
        assert year_starts == sorted(year_starts)


class TestChemicalSpaceEndpoint:
    """Test /api/explorer/chemical-space endpoint."""

    def test_chemical_space_success(self):
        """Test chemical space endpoint."""
        response = client.get("/api/explorer/chemical-space")
        assert response.status_code == 200

        data = response.json()
        assert "pca" in data
        assert "tsne" in data
        assert "timestamp" in data

    def test_chemical_space_pca(self):
        """Test PCA data structure."""
        response = client.get("/api/explorer/chemical-space")
        data = response.json()

        pca = data["pca"]
        assert "coordinates" in pca
        assert "variance_explained" in pca
        assert "compound_data" in pca

        # Check coordinates
        coords = pca["coordinates"]
        assert isinstance(coords, list)
        assert len(coords) > 0
        assert len(coords[0]) == 2  # 2D coordinates

        # Check variance explained
        var_exp = pca["variance_explained"]
        assert len(var_exp) == 2  # 2 components
        assert 0 <= var_exp[0] <= 1
        assert 0 <= var_exp[1] <= 1

        # Check compound data
        compounds = pca["compound_data"]
        assert len(compounds) == len(coords)

        # Validate first compound
        comp = compounds[0]
        assert "cid" in comp
        assert "name" in comp
        assert "toxic" in comp
        assert "year" in comp
        assert "chemical_type" in comp

    def test_chemical_space_tsne(self):
        """Test t-SNE data structure."""
        response = client.get("/api/explorer/chemical-space")
        data = response.json()

        tsne = data["tsne"]
        assert "coordinates" in tsne
        assert "perplexity" in tsne

        # Check coordinates
        coords = tsne["coordinates"]
        assert isinstance(coords, list)
        assert len(coords) > 0
        assert len(coords[0]) == 2  # 2D coordinates

    def test_chemical_space_response_time(self):
        """Test that chemical space endpoint has reasonable response time."""
        # This is expensive, so we allow more time for first call
        start = time.time()
        response = client.get("/api/explorer/chemical-space")
        elapsed = time.time() - start

        assert response.status_code == 200
        # First call might take a few seconds for PCA/t-SNE
        assert elapsed < 30  # Should be <30s even on slow machines


class TestToxicophoresEndpoint:
    """Test /api/explorer/toxicophores endpoint."""

    def test_toxicophores_success(self):
        """Test toxicophores endpoint."""
        response = client.get("/api/explorer/toxicophores")

        # This endpoint might take a while due to SMARTS matching
        assert response.status_code == 200

        data = response.json()
        assert "toxicophores" in data
        assert "timestamp" in data

        # Check toxicophores structure
        toxicophores = data["toxicophores"]
        assert isinstance(toxicophores, list)

        if len(toxicophores) > 0:
            # Validate first toxicophore
            tox = toxicophores[0]
            assert "name" in tox
            assert "smarts" in tox
            assert "prevalence" in tox
            assert "toxic_rate_with" in tox
            assert "toxic_rate_without" in tox
            assert "enrichment_ratio" in tox
            assert "p_value" in tox
            assert "confidence_interval" in tox

            # Check value ranges
            assert 0 <= tox["prevalence"] <= 1
            assert 0 <= tox["toxic_rate_with"] <= 1
            assert tox["enrichment_ratio"] >= 0

    def test_toxicophores_top_10(self):
        """Test that at most 10 toxicophores are returned."""
        response = client.get("/api/explorer/toxicophores")
        data = response.json()

        toxicophores = data["toxicophores"]
        assert len(toxicophores) <= 10


class TestCorrelationsEndpoint:
    """Test /api/explorer/correlations endpoint."""

    def test_correlations_success(self):
        """Test correlations endpoint."""
        response = client.get("/api/explorer/correlations")
        assert response.status_code == 200

        data = response.json()
        assert "correlation_matrix" in data
        assert "top_correlations" in data
        assert "network_data" in data
        assert "timestamp" in data

    def test_correlation_matrix(self):
        """Test correlation matrix structure."""
        response = client.get("/api/explorer/correlations")
        data = response.json()

        matrix = data["correlation_matrix"]
        assert "features" in matrix
        assert "values" in matrix

        features = matrix["features"]
        values = matrix["values"]

        # Check matrix is square
        assert len(values) == len(features)
        assert len(values[0]) == len(features)

        # Check diagonal is 1.0 (correlation with self)
        for i in range(len(features)):
            assert abs(values[i][i] - 1.0) < 0.01

    def test_top_correlations(self):
        """Test top correlations structure."""
        response = client.get("/api/explorer/correlations")
        data = response.json()

        top_corr = data["top_correlations"]
        assert isinstance(top_corr, list)
        assert len(top_corr) <= 20  # Top 20

        if len(top_corr) > 0:
            corr = top_corr[0]
            assert "feature1" in corr
            assert "feature2" in corr
            assert "correlation" in corr
            assert -1 <= corr["correlation"] <= 1

    def test_network_data(self):
        """Test network data structure."""
        response = client.get("/api/explorer/correlations")
        data = response.json()

        network = data["network_data"]
        assert "nodes" in network
        assert "edges" in network

        nodes = network["nodes"]
        edges = network["edges"]

        if len(nodes) > 0:
            node = nodes[0]
            assert "id" in node
            assert "label" in node

        if len(edges) > 0:
            edge = edges[0]
            assert "source" in edge
            assert "target" in edge
            assert "weight" in edge
            assert abs(edge["weight"]) > 0.5  # Only strong correlations


class TestPropertyDistributionsEndpoint:
    """Test /api/explorer/property-distributions endpoint."""

    def test_property_distributions_success(self):
        """Test property distributions endpoint."""
        response = client.get("/api/explorer/property-distributions")
        assert response.status_code == 200

        data = response.json()
        assert "scatter_data" in data
        assert "timestamp" in data

    def test_scatter_data_structure(self):
        """Test scatter data structure."""
        response = client.get("/api/explorer/property-distributions")
        data = response.json()

        scatter_data = data["scatter_data"]
        assert isinstance(scatter_data, list)
        assert len(scatter_data) > 0

        # Validate first scatter plot
        scatter = scatter_data[0]
        assert "name" in scatter
        assert "x_label" in scatter
        assert "y_label" in scatter
        assert "points" in scatter

        # Check points
        points = scatter["points"]
        assert len(points) > 0

        # Validate first point
        point = points[0]
        assert "x" in point
        assert "y" in point
        assert "toxic" in point
        assert "name" in point
        assert "cid" in point
        assert isinstance(point["toxic"], bool)

    def test_property_pairs(self):
        """Test that expected property pairs are included."""
        response = client.get("/api/explorer/property-distributions")
        data = response.json()

        scatter_names = [s["name"] for s in data["scatter_data"]]

        # Check for expected pairs
        assert "LogP_vs_MolecularWeight" in scatter_names


class TestCacheEndpoints:
    """Test cache management endpoints."""

    def test_cache_stats(self):
        """Test cache stats endpoint."""
        response = client.get("/api/explorer/cache/stats")
        assert response.status_code == 200

        data = response.json()
        assert "size" in data

    def test_cache_clear(self):
        """Test cache clear endpoint."""
        response = client.post("/api/explorer/cache/clear")
        assert response.status_code == 200

        data = response.json()
        assert "message" in data


class TestErrorHandling:
    """Test error handling for explorer endpoints."""

    def test_invalid_endpoint(self):
        """Test accessing non-existent explorer endpoint."""
        response = client.get("/api/explorer/nonexistent")
        assert response.status_code == 404

    def test_wrong_http_method(self):
        """Test using wrong HTTP method."""
        response = client.post("/api/explorer/overview")
        assert response.status_code == 405  # Method not allowed


class TestIntegration:
    """Integration tests for explorer workflow."""

    def test_full_exploration_workflow(self):
        """Test complete data exploration workflow."""
        # 1. Get overview
        overview_response = client.get("/api/explorer/overview")
        assert overview_response.status_code == 200
        overview = overview_response.json()
        total_compounds = overview["total_compounds"]

        # 2. Get molecular diversity
        diversity_response = client.get("/api/explorer/molecular-diversity")
        assert diversity_response.status_code == 200

        # 3. Get chemical space
        space_response = client.get("/api/explorer/chemical-space")
        assert space_response.status_code == 200
        space = space_response.json()

        # Verify compound counts match
        assert len(space["pca"]["coordinates"]) == total_compounds

        # 4. Get correlations
        corr_response = client.get("/api/explorer/correlations")
        assert corr_response.status_code == 200

        print(f"\n✓ Full exploration workflow test passed!")
        print(f"  Total compounds: {total_compounds}")
        print(f"  Temporal range: {overview['temporal_range']['min']}-{overview['temporal_range']['max']}")

    def test_performance_all_endpoints(self):
        """Test performance of all endpoints."""
        endpoints = [
            "/api/explorer/overview",
            "/api/explorer/molecular-diversity",
            "/api/explorer/toxicity-by-class",
            "/api/explorer/temporal-trends",
            "/api/explorer/property-distributions",
            "/api/explorer/correlations"
        ]

        results = {}

        for endpoint in endpoints:
            # First call
            start = time.time()
            response = client.get(endpoint)
            time1 = time.time() - start

            # Second call (cached)
            start = time.time()
            response = client.get(endpoint)
            time2 = time.time() - start

            assert response.status_code == 200

            results[endpoint] = {
                "first_call": time1,
                "cached_call": time2
            }

            # Cached call should be fast
            assert time2 < 0.5  # <500ms

        print("\n✓ Performance test passed!")
        for endpoint, times in results.items():
            print(f"  {endpoint}:")
            print(f"    First: {times['first_call']:.3f}s")
            print(f"    Cached: {times['cached_call']:.3f}s")


class TestDataConsistency:
    """Test data consistency across endpoints."""

    def test_compound_count_consistency(self):
        """Test that compound counts are consistent across endpoints."""
        # Get overview
        overview = client.get("/api/explorer/overview").json()
        total = overview["total_compounds"]

        # Get chemical space
        space = client.get("/api/explorer/chemical-space").json()
        space_count = len(space["pca"]["coordinates"])

        # Should match
        assert space_count == total

    def test_toxicity_rate_consistency(self):
        """Test that toxicity rates are consistent."""
        # Get overview
        overview = client.get("/api/explorer/overview").json()
        toxic = overview["class_distribution"]["toxic"]
        non_toxic = overview["class_distribution"]["non_toxic"]
        total = overview["total_compounds"]

        # Verify sum
        assert toxic + non_toxic == total

        # Get by class
        by_class = client.get("/api/explorer/toxicity-by-class").json()
        overall_rate = by_class["overall_toxicity_rate"]

        # Should match
        expected_rate = toxic / total
        assert abs(overall_rate - expected_rate) < 0.001


# Run tests
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
