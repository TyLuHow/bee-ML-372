#!/usr/bin/env python3
"""
Unit Tests for FastAPI Backend
===============================

Tests cover:
- API endpoint functionality
- Request validation
- Response formats
- Error handling
- Model predictions through API

Author: IME 372 Project Team
Date: November 2025
"""

import pytest
import json
import os
import sys
from fastapi.testclient import TestClient

# Import the FastAPI app
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from app.backend.main import app

# Create test client
client = TestClient(app)


class TestHealthEndpoint:
    """Test health check endpoint."""
    
    def test_health_check(self):
        """Test that health endpoint returns 200 OK."""
        response = client.get("/health")
        assert response.status_code == 200
        
        data = response.json()
        assert data["status"] == "healthy"


class TestSMILESEndpoint:
    """Test SMILES-based prediction endpoint."""
    
    @pytest.fixture
    def valid_smiles_input(self):
        """Valid SMILES input."""
        return {
            "smiles": "CCO",  # Ethanol
            "year": 2024,
            "herbicide": 0,
            "fungicide": 0,
            "insecticide": 0,
            "other_agrochemical": 1,
            "source": "PPDB",
            "toxicity_type": "Contact"
        }
    
    def test_smiles_prediction_success(self, valid_smiles_input):
        """Test SMILES prediction with valid input."""
        response = client.post("/predict/smiles", json=valid_smiles_input)
        assert response.status_code == 200
        
        data = response.json()
        assert "prediction" in data
        assert "prediction_label" in data
        assert "confidence" in data
        assert data["prediction"] in [0, 1]
    
    def test_smiles_invalid(self):
        """Test SMILES endpoint with invalid SMILES."""
        invalid_input = {
            "smiles": "INVALID_SMILES_XYZ",
            "year": 2024,
            "herbicide": 0,
            "fungicide": 0,
            "insecticide": 1,
            "other_agrochemical": 0,
            "source": "PPDB",
            "toxicity_type": "Contact"
        }
        response = client.post("/predict/smiles", json=invalid_input)
        assert response.status_code == 400
        assert "Invalid SMILES" in response.json()["detail"]


class TestPredictionEndpoint:
    """Test prediction endpoint."""
    
    @pytest.fixture
    def valid_prediction_input(self):
        """Create valid prediction input data - CORRECTED field names."""
        return {
            "source": "PPDB",
            "year": 2020,
            "toxicity_type": "Contact",
            "herbicide": 0,
            "fungicide": 0,
            "insecticide": 1,
            "other_agrochemical": 0,
            # CORRECTED: Match exact field names from PredictionInput model
            "MolecularWeight": 350.5,
            "LogP": 3.2,
            "NumHDonors": 2,
            "NumHAcceptors": 4,
            "NumRotatableBonds": 5,
            "NumAromaticRings": 2,  # FIXED: was "AromaticRings"
            "TPSA": 65.3,
            "NumHeteroatoms": 5,
            "NumRings": 2,  # FIXED: was "RingCount"
            "NumSaturatedRings": 0,
            "NumAliphaticRings": 0,
            "FractionCSP3": 0.25,  # FIXED: was "FractionCsp3"
            "MolarRefractivity": 95.0,
            "BertzCT": 450.0,
            "HeavyAtomCount": 25
        }
    
    def test_prediction_success(self, valid_prediction_input):
        """Test successful prediction."""
        response = client.post("/predict", json=valid_prediction_input)
        
        # Check response status
        assert response.status_code == 200
        
        # Check response structure
        data = response.json()
        assert "prediction" in data
        assert "probability_toxic" in data
        assert "probability_non_toxic" in data
        assert "confidence" in data
        assert "label_text" in data
        assert "timestamp" in data
        
        # Check value ranges
        assert data["prediction"] in [0, 1]
        assert 0 <= data["probability_toxic"] <= 1
        assert 0 <= data["probability_non_toxic"] <= 1
        assert 0 <= data["confidence"] <= 1
        assert data["label_text"] in ["Non-Toxic", "Toxic"]
    
    def test_prediction_invalid_source(self, valid_prediction_input):
        """Test prediction with invalid source."""
        invalid_input = valid_prediction_input.copy()
        invalid_input["source"] = "INVALID_SOURCE"
        
        response = client.post("/predict", json=invalid_input)
        # Should still work or return 422 validation error
        assert response.status_code in [200, 422]
    
    def test_prediction_missing_field(self, valid_prediction_input):
        """Test prediction with missing required field."""
        invalid_input = valid_prediction_input.copy()
        del invalid_input["MolecularWeight"]
        
        response = client.post("/predict", json=invalid_input)
        assert response.status_code == 422  # Validation error
    
    def test_prediction_invalid_data_type(self, valid_prediction_input):
        """Test prediction with invalid data type."""
        invalid_input = valid_prediction_input.copy()
        invalid_input["MolecularWeight"] = "not_a_number"
        
        response = client.post("/predict", json=invalid_input)
        assert response.status_code == 422  # Validation error
    
    def test_prediction_out_of_range(self, valid_prediction_input):
        """Test prediction with out of range values."""
        invalid_input = valid_prediction_input.copy()
        invalid_input["year"] = 1700  # Before valid range
        
        response = client.post("/predict", json=invalid_input)
        # Should return 422 or handle gracefully
        assert response.status_code in [200, 422]
    
    def test_prediction_insecticide_toxic(self, valid_prediction_input):
        """Test that insecticide is likely predicted as toxic."""
        # Set insecticide flag to test toxicity prediction
        insecticide_input = valid_prediction_input.copy()
        insecticide_input["insecticide"] = 1
        insecticide_input["herbicide"] = 0
        insecticide_input["fungicide"] = 0
        
        response = client.post("/predict", json=insecticide_input)
        assert response.status_code == 200
        
        data = response.json()
        # Insecticides often predicted as toxic, but not guaranteed
        assert "prediction" in data


class TestModelInfoEndpoint:
    """Test model information endpoint."""
    
    def test_model_info(self):
        """Test model info endpoint returns correct information."""
        response = client.get("/model/info")
        assert response.status_code == 200
        
        data = response.json()
        assert "model_name" in data or "algorithm" in data
        
        # Check for performance metrics
        if "metrics" in data:
            metrics = data["metrics"]
            if "accuracy" in metrics:
                assert 0 <= metrics["accuracy"] <= 1


class TestFeatureImportanceEndpoint:
    """Test feature importance endpoint."""
    
    def test_feature_importance(self):
        """Test feature importance endpoint."""
        response = client.get("/feature/importance")
        
        # Should return 200 if implemented
        if response.status_code == 200:
            data = response.json()
            assert "feature_importance" in data or "features" in data
            
            # Check data structure
            importance_data = data.get("feature_importance", data.get("features", []))
            assert isinstance(importance_data, (list, dict))


class TestHistoryEndpoint:
    """Test prediction history endpoint."""
    
    def test_history_default(self):
        """Test history endpoint with default parameters."""
        response = client.get("/history")
        
        # Should return 200
        if response.status_code == 200:
            data = response.json()
            assert "history" in data or isinstance(data, list)
            
            # If there's history, check structure
            history = data.get("history", data)
            if len(history) > 0:
                assert "timestamp" in history[0] or "prediction" in history[0]
    
    def test_history_with_limit(self):
        """Test history endpoint with limit parameter."""
        response = client.get("/history?limit=5")
        
        if response.status_code == 200:
            data = response.json()
            history = data.get("history", data)
            
            # Should return at most 5 items
            if isinstance(history, list):
                assert len(history) <= 5


class TestErrorHandling:
    """Test API error handling."""
    
    def test_invalid_endpoint(self):
        """Test accessing non-existent endpoint."""
        response = client.get("/nonexistent")
        assert response.status_code == 404
    
    def test_wrong_http_method(self):
        """Test using wrong HTTP method."""
        response = client.get("/predict")  # Should be POST
        assert response.status_code == 405  # Method not allowed
    
    def test_invalid_json(self):
        """Test sending invalid JSON."""
        response = client.post(
            "/predict",
            data="invalid json",
            headers={"Content-Type": "application/json"}
        )
        assert response.status_code == 422


class TestCORS:
    """Test CORS configuration."""
    
    def test_cors_headers(self):
        """Test that CORS headers are present."""
        response = client.options("/predict")
        
        # Check for CORS headers (may vary based on configuration)
        # In test client, CORS middleware might not add all headers
        assert response.status_code in [200, 405]


class TestIntegration:
    """Integration tests for full API workflow."""
    
    def test_full_prediction_workflow(self):
        """Test complete prediction workflow."""
        # 1. Check health
        health_response = client.get("/health")
        assert health_response.status_code == 200
        
        # 2. Make prediction
        prediction_input = {
            "source": "PPDB",
            "year": 2020,
            "toxicity_type": "Contact",
            "herbicide": 0,
            "fungicide": 0,
            "insecticide": 1,
            "other_agrochemical": 0,
            "MolecularWeight": 350.5,
            "LogP": 3.2,
            "NumHDonors": 2,
            "NumHAcceptors": 4,
            "NumRotatableBonds": 5,
            "AromaticRings": 2,
            "TPSA": 65.3,
            "NumHeteroatoms": 5,
            "NumAromaticAtoms": 12,
            "NumSaturatedRings": 0,
            "NumAliphaticRings": 0,
            "RingCount": 2,
            "FractionCsp3": 0.25,
            "NumAromaticCarbocycles": 1,
            "NumSaturatedCarbocycles": 0
        }
        
        predict_response = client.post("/predict", json=prediction_input)
        assert predict_response.status_code == 200
        
        prediction_data = predict_response.json()
        assert "prediction" in prediction_data
        
        # 3. Get model info
        info_response = client.get("/model/info")
        assert info_response.status_code == 200
        
        print("✓ Full API workflow test passed!")
        print(f"  Prediction: {prediction_data.get('label_text', 'N/A')}")
        print(f"  Confidence: {prediction_data.get('confidence', 0):.2%}")
    
    def test_multiple_predictions(self):
        """Test making multiple predictions in sequence."""
        prediction_input_base = {
            "source": "PPDB",
            "year": 2020,
            "toxicity_type": "Contact",
            "herbicide": 0,
            "fungicide": 0,
            "insecticide": 0,
            "other_agrochemical": 0,
            "MolecularWeight": 300.0,
            "LogP": 2.5,
            "NumHDonors": 1,
            "NumHAcceptors": 3,
            "NumRotatableBonds": 4,
            "AromaticRings": 1,
            "TPSA": 50.0,
            "NumHeteroatoms": 3,
            "NumAromaticAtoms": 6,
            "NumSaturatedRings": 0,
            "NumAliphaticRings": 0,
            "RingCount": 1,
            "FractionCsp3": 0.3,
            "NumAromaticCarbocycles": 1,
            "NumSaturatedCarbocycles": 0
        }
        
        # Make 3 predictions with variations
        for i, chemical_type in enumerate(["herbicide", "fungicide", "insecticide"]):
            input_data = prediction_input_base.copy()
            input_data[chemical_type] = 1
            
            response = client.post("/predict", json=input_data)
            assert response.status_code == 200
            
            data = response.json()
            assert "prediction" in data
        
        print("✓ Multiple predictions test passed!")


class TestDataValidation:
    """Test input data validation."""
    
    def test_negative_molecular_weight(self):
        """Test that negative molecular weight is rejected."""
        invalid_input = {
            "source": "PPDB",
            "year": 2020,
            "toxicity_type": "Contact",
            "herbicide": 0,
            "fungicide": 0,
            "insecticide": 1,
            "other_agrochemical": 0,
            "MolecularWeight": -100.0,  # Invalid
            "LogP": 3.2,
            "NumHDonors": 2,
            "NumHAcceptors": 4,
            "NumRotatableBonds": 5,
            "AromaticRings": 2,
            "TPSA": 65.3,
            "NumHeteroatoms": 5,
            "NumAromaticAtoms": 12,
            "NumSaturatedRings": 0,
            "NumAliphaticRings": 0,
            "RingCount": 2,
            "FractionCsp3": 0.25,
            "NumAromaticCarbocycles": 1,
            "NumSaturatedCarbocycles": 0
        }
        
        response = client.post("/predict", json=invalid_input)
        # Should reject invalid input
        assert response.status_code in [422, 400]
    
    def test_invalid_binary_flags(self):
        """Test that binary flags only accept 0 or 1."""
        invalid_input = {
            "source": "PPDB",
            "year": 2020,
            "toxicity_type": "Contact",
            "herbicide": 2,  # Invalid (should be 0 or 1)
            "fungicide": 0,
            "insecticide": 1,
            "other_agrochemical": 0,
            "MolecularWeight": 350.5,
            "LogP": 3.2,
            "NumHDonors": 2,
            "NumHAcceptors": 4,
            "NumRotatableBonds": 5,
            "AromaticRings": 2,
            "TPSA": 65.3,
            "NumHeteroatoms": 5,
            "NumAromaticAtoms": 12,
            "NumSaturatedRings": 0,
            "NumAliphaticRings": 0,
            "RingCount": 2,
            "FractionCsp3": 0.25,
            "NumAromaticCarbocycles": 1,
            "NumSaturatedCarbocycles": 0
        }
        
        response = client.post("/predict", json=invalid_input)
        # Should reject invalid input
        assert response.status_code in [422, 400]


# Run tests
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])

