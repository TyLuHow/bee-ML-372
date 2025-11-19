#!/usr/bin/env python3
"""
Molecular Feature Extraction Module
====================================

Generate molecular descriptors from SMILES using RDKit.
Supports conversion of SMILES strings to standardized molecular descriptors
for toxicity prediction.

Author: IME 372 Project Team
Date: November 2025
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import pandas as pd
import numpy as np
from typing import Optional, Dict, List
import warnings
warnings.filterwarnings('ignore')


class MolecularFeaturizer:
    """
    Generate molecular descriptors from SMILES using RDKit.
    
    Features:
    - SMILES parsing and validation
    - Molecular descriptor calculation (15 descriptors)
    - Batch processing of multiple SMILES
    - Error handling for invalid structures
    """
    
    def __init__(self):
        """Initialize featurizer with descriptor definitions."""
        self.descriptor_functions = {
            'MolecularWeight': Descriptors.MolWt,
            'LogP': Crippen.MolLogP,
            'NumHDonors': Lipinski.NumHDonors,
            'NumHAcceptors': Lipinski.NumHAcceptors,
            'NumRotatableBonds': Lipinski.NumRotatableBonds,
            'NumAromaticRings': Lipinski.NumAromaticRings,
            'TPSA': rdMolDescriptors.CalcTPSA,
            'NumHeteroatoms': Lipinski.NumHeteroatoms,
            'NumRings': Lipinski.RingCount,
            'NumSaturatedRings': Lipinski.NumSaturatedRings,
            'NumAliphaticRings': Lipinski.NumAliphaticRings,
            'FractionCSP3': Lipinski.FractionCSP3,
            'MolarRefractivity': Crippen.MolMR,
            'BertzCT': Descriptors.BertzCT,
            'HeavyAtomCount': Lipinski.HeavyAtomCount
        }
        
        self.feature_names = list(self.descriptor_functions.keys())
    
    def smiles_to_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Parse SMILES to RDKit molecule object.
        
        Args:
            smiles: SMILES string
            
        Returns:
            RDKit Mol object or None if invalid
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"X Invalid SMILES: {smiles}")
                return None
            return mol
        except Exception as e:
            print(f"X SMILES parsing error for '{smiles}': {e}")
            return None
    
    def calculate_descriptors(self, mol: Chem.Mol) -> Dict[str, float]:
        """
        Calculate all molecular descriptors for a molecule.
        
        Args:
            mol: RDKit Mol object
            
        Returns:
            Dictionary of descriptor name: value pairs
        """
        descriptors = {}
        for name, func in self.descriptor_functions.items():
            try:
                value = func(mol)
                # Handle potential None/NaN values
                if value is None or (isinstance(value, float) and np.isnan(value)):
                    descriptors[name] = 0.0
                else:
                    descriptors[name] = float(value)
            except Exception as e:
                print(f"X Error calculating {name}: {e}")
                descriptors[name] = 0.0
        
        return descriptors
    
    def smiles_to_descriptors(self, smiles: str) -> Optional[Dict[str, float]]:
        """
        Convert SMILES to descriptor dictionary (one-step convenience).
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary of descriptors or None if invalid SMILES
        """
        mol = self.smiles_to_mol(smiles)
        if mol is None:
            return None
        return self.calculate_descriptors(mol)
    
    def batch_smiles_to_dataframe(
        self, 
        smiles_list: List[str],
        include_invalid: bool = False
    ) -> pd.DataFrame:
        """
        Convert multiple SMILES to descriptor DataFrame.
        
        Args:
            smiles_list: List of SMILES strings
            include_invalid: If True, include rows for invalid SMILES with NaN values
            
        Returns:
            DataFrame with descriptors (rows=molecules, columns=descriptors)
        """
        descriptors_list = []
        
        for i, smiles in enumerate(smiles_list):
            desc = self.smiles_to_descriptors(smiles)
            
            if desc is not None:
                desc['SMILES'] = smiles
                desc['index'] = i
                descriptors_list.append(desc)
            elif include_invalid:
                # Add row with NaN for invalid SMILES
                invalid_desc = {name: np.nan for name in self.feature_names}
                invalid_desc['SMILES'] = smiles
                invalid_desc['index'] = i
                descriptors_list.append(invalid_desc)
        
        if not descriptors_list:
            return pd.DataFrame(columns=self.feature_names + ['SMILES', 'index'])
        
        df = pd.DataFrame(descriptors_list)
        
        print(f"OK Processed {len(descriptors_list)}/{len(smiles_list)} SMILES")
        if len(descriptors_list) < len(smiles_list):
            print(f"  WARNING: {len(smiles_list) - len(descriptors_list)} invalid SMILES skipped")
        
        return df
    
    def validate_smiles(self, smiles: str) -> Dict[str, any]:
        """
        Validate SMILES and return detailed information.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary with validation results
        """
        mol = self.smiles_to_mol(smiles)
        
        if mol is None:
            return {
                'valid': False,
                'smiles': smiles,
                'error': 'Invalid SMILES format',
                'canonical_smiles': None,
                'num_atoms': None,
                'num_bonds': None
            }
        
        try:
            canonical_smiles = Chem.MolToSmiles(mol)
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            
            return {
                'valid': True,
                'smiles': smiles,
                'canonical_smiles': canonical_smiles,
                'num_atoms': num_atoms,
                'num_bonds': num_bonds,
                'error': None
            }
        except Exception as e:
            return {
                'valid': False,
                'smiles': smiles,
                'error': str(e),
                'canonical_smiles': None,
                'num_atoms': None,
                'num_bonds': None
            }
    
    def get_feature_info(self) -> pd.DataFrame:
        """
        Get information about available molecular descriptors.
        
        Returns:
            DataFrame with feature names and descriptions
        """
        feature_info = [
            {'name': 'MolecularWeight', 'description': 'Molecular weight in Daltons', 'type': 'continuous'},
            {'name': 'LogP', 'description': 'Octanol-water partition coefficient (lipophilicity)', 'type': 'continuous'},
            {'name': 'NumHDonors', 'description': 'Number of hydrogen bond donors', 'type': 'integer'},
            {'name': 'NumHAcceptors', 'description': 'Number of hydrogen bond acceptors', 'type': 'integer'},
            {'name': 'NumRotatableBonds', 'description': 'Number of rotatable bonds (flexibility)', 'type': 'integer'},
            {'name': 'NumAromaticRings', 'description': 'Number of aromatic rings', 'type': 'integer'},
            {'name': 'TPSA', 'description': 'Topological polar surface area', 'type': 'continuous'},
            {'name': 'NumHeteroatoms', 'description': 'Number of non-carbon atoms', 'type': 'integer'},
            {'name': 'NumRings', 'description': 'Total number of rings', 'type': 'integer'},
            {'name': 'NumSaturatedRings', 'description': 'Number of saturated rings', 'type': 'integer'},
            {'name': 'NumAliphaticRings', 'description': 'Number of aliphatic rings', 'type': 'integer'},
            {'name': 'FractionCSP3', 'description': 'Fraction of sp3 hybridized carbons', 'type': 'continuous'},
            {'name': 'MolarRefractivity', 'description': 'Molar refractivity', 'type': 'continuous'},
            {'name': 'BertzCT', 'description': 'Bertz molecular complexity index', 'type': 'continuous'},
            {'name': 'HeavyAtomCount', 'description': 'Number of heavy (non-hydrogen) atoms', 'type': 'integer'}
        ]
        
        return pd.DataFrame(feature_info)


# Example usage
if __name__ == "__main__":
    print("="*80)
    print("MOLECULAR FEATURIZER - Example Usage")
    print("="*80)
    
    # Initialize featurizer
    featurizer = MolecularFeaturizer()
    
    # Example SMILES
    example_smiles = [
        "CCO",  # Ethanol
        "CC(C)C",  # Isobutane
        "c1ccccc1",  # Benzene
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "INVALID_SMILES"  # Invalid
    ]
    
    print("\n1. Single SMILES processing:")
    for smiles in example_smiles[:3]:
        descriptors = featurizer.smiles_to_descriptors(smiles)
        if descriptors:
            print(f"\n  SMILES: {smiles}")
            print(f"  MolWt: {descriptors['MolecularWeight']:.2f}")
            print(f"  LogP: {descriptors['LogP']:.2f}")
            print(f"  NumHDonors: {descriptors['NumHDonors']}")
    
    print("\n2. Batch processing:")
    df = featurizer.batch_smiles_to_dataframe(example_smiles)
    print(f"\n  Result shape: {df.shape}")
    print(f"\n  Columns: {list(df.columns)}")
    print(f"\n  First 3 rows:")
    print(df.head(3))
    
    print("\n3. Feature information:")
    info_df = featurizer.get_feature_info()
    print(info_df)
    
    print("\nOK Example completed successfully!")


