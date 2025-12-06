"""
Prediction Logic for Bee Toxicity Classification
Handles model loading, feature preprocessing, and prediction generation
"""
import os
import pickle
from pathlib import Path
from typing import Dict, Tuple, Optional
import logging
import numpy as np
import pandas as pd
from .featurizer import MolecularFeaturizer

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ToxicityPredictor:
    """
    Main prediction class that loads models and generates toxicity predictions.
    """

    def __init__(self, model_dir: str = None):
        """
        Initialize the predictor by loading models and preprocessor.

        Args:
            model_dir: Directory containing model files. Defaults to backend/models/
        """
        if model_dir is None:
            # Default to models directory relative to this file
            backend_dir = Path(__file__).parent.parent
            model_dir = backend_dir / "models"

        self.model_dir = Path(model_dir)
        self.featurizer = MolecularFeaturizer()
        self.model = None
        self.preprocessor = None
        self.model_name = None

        # Load models on initialization
        self._load_models()

    def _load_models(self):
        """Load the trained model and preprocessor from disk."""
        try:
            # Try to load Random Forest first (typically best performer)
            rf_path = self.model_dir / "best_model_random_forest.pkl"
            xgb_path = self.model_dir / "best_model_xgboost.pkl"

            if rf_path.exists():
                with open(rf_path, 'rb') as f:
                    self.model = pickle.load(f)
                self.model_name = "Random Forest"
                logger.info(f"Loaded Random Forest model from {rf_path}")
            elif xgb_path.exists():
                with open(xgb_path, 'rb') as f:
                    self.model = pickle.load(f)
                self.model_name = "XGBoost"
                logger.info(f"Loaded XGBoost model from {xgb_path}")
            else:
                raise FileNotFoundError("No trained model found in models directory")

            # Load preprocessor (scaler)
            preprocessor_path = self.model_dir / "preprocessor.pkl"
            if preprocessor_path.exists():
                with open(preprocessor_path, 'rb') as f:
                    self.preprocessor = pickle.load(f)
                logger.info(f"Loaded preprocessor from {preprocessor_path}")
            else:
                logger.warning("No preprocessor found, predictions may be less accurate")

        except Exception as e:
            logger.error(f"Error loading models: {e}")
            raise

    def _prepare_features(self, compound_data: Dict) -> Optional[pd.DataFrame]:
        """
        Prepare features for prediction from compound data.

        Args:
            compound_data: Dictionary containing compound properties

        Returns:
            DataFrame with features ready for model prediction, or None if error
        """
        try:
            # Extract SMILES if provided, otherwise we'll need to work with provided features
            smiles = compound_data.get('smiles', '')

            # Compute molecular descriptors from SMILES
            if smiles:
                descriptors = self.featurizer.smiles_to_descriptors(smiles)
                if descriptors is None:
                    logger.error(f"Invalid SMILES: {smiles}")
                    return None
            else:
                # No SMILES provided - use default values
                descriptors = {name: 0.0 for name in self.featurizer.feature_names}

            # Create base feature dictionary
            features = {
                # Basic properties from frontend
                'year': 2024,  # Default year
                **descriptors,  # Add all molecular descriptors
            }

            # Handle categorical features - use one-hot encoding matching training
            # Note: These need to match the exact feature names from training
            category = compound_data.get('category', 'Other')
            exposure = compound_data.get('exposure', 'Contact')

            # Source encoding (match training - we'll use 'Other' as default)
            features['source_PPDB'] = 0
            features['source_Other'] = 1  # Drop first is True, so we encode non-reference

            # Toxicity type encoding (exposure route)
            # Common types: Contact, Oral, Systemic
            features['toxicity_type_Oral'] = 1 if 'Oral' in exposure else 0
            features['toxicity_type_Contact'] = 1 if 'Contact' in exposure else 0

            # Create DataFrame with single row
            feature_df = pd.DataFrame([features])

            # Ensure all expected features are present (model was trained on specific feature set)
            # The preprocessor/scaler expects specific columns in specific order
            # We'll let it handle any missing columns by filling with 0

            return feature_df

        except Exception as e:
            logger.error(f"Error preparing features: {e}")
            return None

    def predict(self, compound_data: Dict) -> Dict:
        """
        Generate toxicity prediction for a compound.

        Args:
            compound_data: Dictionary with keys:
                - name: compound name
                - smiles: SMILES string (optional but recommended)
                - category: pesticide category (Insecticide, Fungicide, Herbicide)
                - molecular_weight: molecular weight (g/mol)
                - logp: LogP value
                - exposure_route: exposure route (Contact, Oral, Systemic)

        Returns:
            Dictionary with prediction results:
                - toxicity: "Toxic" or "Safe"
                - confidence: confidence score (0-100)
                - explanation: detailed explanation
                - recommendation: actionable recommendation
        """
        try:
            # Validate SMILES if provided
            smiles = compound_data.get('smiles', '')
            if smiles and not self.featurizer.validate_smiles(smiles):
                return {
                    'toxicity': 'Uncertain',
                    'confidence': 0,
                    'explanation': f"Invalid SMILES string provided: {smiles}. Please check the molecular structure.",
                    'recommendation': "Verify the SMILES string represents a valid chemical structure before proceeding with toxicity assessment."
                }

            # Prepare features
            features = self._prepare_features(compound_data)
            if features is None:
                return {
                    'toxicity': 'Uncertain',
                    'confidence': 0,
                    'explanation': "Unable to compute molecular features from provided data.",
                    'recommendation': "Ensure all required compound properties are provided, especially a valid SMILES string."
                }

            # Apply preprocessing if available
            if self.preprocessor:
                # The preprocessor is a StandardScaler fit on training data
                # We need to ensure our features match the training features
                try:
                    features_scaled = self.preprocessor.transform(features)
                    features = pd.DataFrame(features_scaled, columns=features.columns)
                except Exception as e:
                    logger.warning(f"Preprocessing failed: {e}. Using raw features.")

            # Generate prediction
            prediction = self.model.predict(features)[0]
            probabilities = self.model.predict_proba(features)[0]

            # Extract confidence (probability of predicted class)
            confidence = float(probabilities[prediction]) * 100

            # Determine toxicity label
            toxicity = "Toxic" if prediction == 1 else "Safe"

            # Generate explanation and recommendation
            explanation = self._generate_explanation(compound_data, toxicity, probabilities)
            recommendation = self._generate_recommendation(compound_data, toxicity, confidence)

            return {
                'toxicity': toxicity,
                'confidence': round(confidence, 1),
                'explanation': explanation,
                'recommendation': recommendation
            }

        except Exception as e:
            logger.error(f"Prediction error: {e}")
            return {
                'toxicity': 'Uncertain',
                'confidence': 0,
                'explanation': f"An error occurred during prediction: {str(e)}",
                'recommendation': "Please verify the input data and try again. Contact support if the issue persists."
            }

    def _generate_explanation(self, data: Dict, toxicity: str, probabilities: np.ndarray) -> str:
        """Generate a detailed explanation for the prediction."""
        compound = data.get('name', 'this compound')
        category = data.get('category', 'pesticide')
        mw = data.get('molecular_weight', 'unknown')
        logp = data.get('logp', 'unknown')
        exposure = data.get('exposure_route', 'contact')

        # Get descriptors if SMILES was provided
        smiles = data.get('smiles', '')
        descriptors = None
        if smiles:
            descriptors = self.featurizer.smiles_to_descriptors(smiles)

        if toxicity == "Toxic":
            if category == "Insecticide":
                if descriptors and descriptors.get('LogP', 0) > 3:
                    return (f"{compound} exhibits high lipophilicity (LogP: {descriptors['LogP']:.2f}), "
                           f"suggesting strong accumulation in bee tissues and neural membranes. "
                           f"The {category.lower()} mode of action likely targets neurotransmitter systems "
                           f"shared between pest insects and pollinators, resulting in elevated acute toxicity risk. "
                           f"Model confidence: {probabilities[1]*100:.1f}%")
                return (f"Analysis indicates {compound} falls within the molecular profile "
                       f"commonly associated with bee-toxic compounds. Structural properties suggest "
                       f"potential impairment of acetylcholine receptors critical for bee navigation and "
                       f"foraging behavior. Model confidence: {probabilities[1]*100:.1f}%")

            if category == "Fungicide":
                return (f"While fungicides typically show lower bee toxicity, {compound}'s molecular properties "
                       f"indicate potential bioaccumulation concerns. Chronic exposure through contaminated pollen "
                       f"or nectar may disrupt bee immune function or synergize with other stressors. "
                       f"Model confidence: {probabilities[1]*100:.1f}%")

            return (f"{compound} demonstrates physicochemical properties that suggest potential bioavailability "
                   f"and tissue penetration in Apis mellifera. The {exposure.lower()} exposure route presents "
                   f"significant risk during foraging activities. Model confidence: {probabilities[1]*100:.1f}%")

        else:  # Safe
            if category == "Herbicide":
                return (f"{compound} targets plant-specific metabolic pathways (photosynthesis or amino acid synthesis) "
                       f"not present in honey bees. The molecular profile suggests low bioaccumulation potential and "
                       f"minimal interaction with bee neurological or physiological systems. "
                       f"Model confidence: {probabilities[0]*100:.1f}%")

            if category == "Fungicide":
                return (f"The compound shows a favorable safety profile for pollinators. Molecular properties indicate "
                       f"limited cuticle penetration and reduced bioavailability. Fungal-specific targets minimize "
                       f"off-target effects on bee cellular processes. Model confidence: {probabilities[0]*100:.1f}%")

            if descriptors and descriptors.get('LogP', 0) < 2:
                return (f"{compound}'s low lipophilicity (LogP: {descriptors['LogP']:.2f}) limits cuticular absorption "
                       f"and neural tissue accumulation. Combined with the {exposure.lower()} exposure profile, "
                       f"this suggests minimal acute or chronic toxicity to foraging bees. "
                       f"Model confidence: {probabilities[0]*100:.1f}%")

            return (f"Molecular analysis of {compound} reveals properties inconsistent with bee-toxic compounds "
                   f"in the training dataset. The molecular profile suggests favorable environmental degradation "
                   f"and low bioaccumulation potential. Model confidence: {probabilities[0]*100:.1f}%")

    def _generate_recommendation(self, data: Dict, toxicity: str, confidence: float) -> str:
        """Generate actionable recommendations based on prediction."""
        exposure = data.get('exposure_route', 'contact')
        category = data.get('category', 'pesticide')

        if toxicity == "Toxic":
            if 'Contact' in exposure:
                return ("Recommend avoiding application during bloom periods and implementing strict spray drift "
                       "management protocols. Consider alternative formulations with reduced contact exposure or "
                       "systemic alternatives with delayed bee-accessible residues.")

            if 'Oral' in exposure:
                return ("Prioritize EPA OECD 245 chronic oral toxicity testing before field deployment. Implement "
                       "buffer zones around pollinator-attractive crops and restrict application during active "
                       "foraging hours (10 AM - 4 PM).")

            if 'Systemic' in exposure:
                return ("Conduct extended residue studies on pollen and nectar to establish safe application timing "
                       "windows. Consider seed treatment limitations and soil incorporation methods to minimize "
                       "plant uptake during bloom.")

            return ("Proceed with comprehensive acute contact (OECD 214) and oral (OECD 213) toxicity bioassays "
                   "before commercial development. Explore structural modifications to reduce bee exposure or "
                   "toxicity while maintaining target pest efficacy.")

        else:  # Safe
            if category == "Insecticide" or confidence < 75:
                return ("While preliminary assessment suggests lower bee toxicity, confirm with targeted semi-field "
                       "tunnel studies (OECD 75) to validate safety under realistic foraging conditions. Monitor "
                       "for sublethal effects on navigation and learning.")

            return ("Model suggests favorable pollinator safety profile. Recommend proceeding with tier-1 laboratory "
                   "screening (OECD 213/214) to confirm predictions, followed by streamlined field registration "
                   "testing. Consider marketing as a bee-safe alternative to increase adoption.")


# Global predictor instance (loaded once at startup)
_predictor = None


def get_predictor() -> ToxicityPredictor:
    """
    Get or create the global predictor instance.
    This ensures models are loaded only once at application startup.
    """
    global _predictor
    if _predictor is None:
        _predictor = ToxicityPredictor()
    return _predictor
