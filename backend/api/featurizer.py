"""
Molecular Feature Engineering using RDKit
Computes molecular descriptors from SMILES strings
"""
from typing import Dict, Optional, List
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors


class MolecularFeaturizer:
    """
    Generate molecular descriptors from SMILES using RDKit.
    This class replicates the feature engineering pipeline from the training notebook.
    """

    def __init__(self):
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

    def smiles_to_descriptors(self, smiles: str) -> Optional[Dict[str, float]]:
        """
        Convert a SMILES string to a dictionary of molecular descriptors.

        Args:
            smiles: SMILES string representation of the molecule

        Returns:
            Dictionary of descriptor name -> value, or None if SMILES is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        descriptors = {}
        for name, func in self.descriptor_functions.items():
            try:
                value = func(mol)
                descriptors[name] = float(value) if value is not None else 0.0
            except Exception:
                descriptors[name] = 0.0

        return descriptors

    def batch_smiles_to_dataframe(self, smiles_list: List[str]) -> pd.DataFrame:
        """
        Convert a list of SMILES to a DataFrame of molecular descriptors.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            DataFrame with molecular descriptors (rows with invalid SMILES are excluded)
        """
        results = []
        for i, smiles in enumerate(smiles_list):
            desc = self.smiles_to_descriptors(smiles)
            if desc:
                results.append(desc)
            if (i + 1) % 200 == 0:
                print(f"  Processed {i + 1}/{len(smiles_list)} SMILES...")

        print(f"Processed {len(results)}/{len(smiles_list)} SMILES successfully")
        return pd.DataFrame(results)

    def validate_smiles(self, smiles: str) -> bool:
        """
        Check if a SMILES string is valid.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
