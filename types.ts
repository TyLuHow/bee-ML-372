// Complete chemical data interface with all molecular descriptors
export interface ChemicalData {
  // Basic Information
  name?: string;
  CID?: number;
  CAS?: string;
  smiles?: string;

  // Data Source & Classification
  source?: 'ECOTOX' | 'PPDB' | 'BPDB' | 'Unknown' | string;
  year?: number;
  toxicity_type?: 'Contact' | 'Oral' | 'Systemic' | 'Other' | string;
  exposure?: string; // Alias for toxicity_type (backwards compatibility)

  // Pesticide Classification (binary flags)
  category?: string; // Human-readable category
  insecticide?: number; // 0 or 1
  herbicide?: number; // 0 or 1
  fungicide?: number; // 0 or 1
  other_agrochemical?: number; // 0 or 1

  // Molecular Descriptors (all 15+ descriptors)
  MolecularWeight?: number;
  mw?: number; // Alias for MolecularWeight (backwards compatibility)
  LogP?: number;
  logP?: number; // Alias for LogP (backwards compatibility)
  NumHDonors?: number;
  NumHAcceptors?: number;
  NumRotatableBonds?: number;
  NumAromaticRings?: number;
  AromaticRings?: number; // Alias
  TPSA?: number;
  NumHeteroatoms?: number;
  NumRings?: number;
  RingCount?: number; // Alias
  NumSaturatedRings?: number;
  NumAliphaticRings?: number;
  FractionCSP3?: number;
  FractionCsp3?: number; // Alias
  MolarRefractivity?: number;
  BertzCT?: number;
  HeavyAtomCount?: number;
  NumAromaticAtoms?: number;
  NumAromaticCarbocycles?: number;
  NumSaturatedCarbocycles?: number;
}

// Backend API prediction result
export interface PredictionResult {
  toxicity?: 'Toxic' | 'Safe' | 'Uncertain';
  prediction?: number; // 0 or 1
  label_text?: string;
  probability_toxic?: number;
  probability_non_toxic?: number;
  confidence: number;
  explanation?: string;
  recommendation?: string;
  timestamp?: string;
}

export type TabView = 'engine' | 'scenarios' | 'science';

// Scenario definition for real-world case studies
export interface Scenario {
  id: string;
  title: string;
  icon: string;
  userType: string;
  role?: string; // Backwards compatibility
  description?: string; // Backwards compatibility
  context: string;
  stakes: string;
  compounds: CompoundData[];
}

// Compound data for scenarios
export interface CompoundData {
  name: string;
  description: string;
  expectedResult: 'TOXIC' | 'NON-TOXIC';
  inputs: Partial<ChemicalData>;
}