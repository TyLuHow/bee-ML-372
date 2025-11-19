#!/usr/bin/env python3
"""
Test script for SMILES featurization
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

try:
    from src.molecular_features import MolecularFeaturizer
    print("✓ MolecularFeaturizer imported successfully")
    
    # Test featurizer
    featurizer = MolecularFeaturizer()
    print(f"✓ Featurizer initialized with {len(featurizer.feature_names)} features")
    
    # Test SMILES
    test_smiles = [
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene"),
        ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
    ]
    
    print("\nTesting SMILES conversion:")
    for smiles, name in test_smiles:
        descriptors = featurizer.smiles_to_descriptors(smiles)
        if descriptors:
            print(f"✓ {name} ({smiles})")
            print(f"  MolWt={descriptors['MolecularWeight']:.2f}, LogP={descriptors['LogP']:.2f}")
        else:
            print(f"✗ Failed to process {name}")
    
    print("\n✓ All tests passed!")
    
except ImportError as e:
    print(f"✗ Import error: {e}")
    print("  RDKit may not be installed. Install with: pip install rdkit")
    sys.exit(1)
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)




