export interface ValidationRule {
  min?: number
  max?: number
  required?: boolean
  message?: string
}

export const validationRules: Record<string, ValidationRule> = {
  year: {
    min: 1800,
    max: 2030,
    required: true,
    message: 'Year must be between 1800 and 2030',
  },
  MolecularWeight: {
    min: 0,
    max: 2000,
    required: true,
    message: 'Molecular weight must be between 0 and 2000',
  },
  LogP: {
    min: -10,
    max: 20,
    required: true,
    message: 'LogP must be between -10 and 20',
  },
  NumHDonors: {
    min: 0,
    max: 20,
    required: true,
    message: 'Number of H-donors must be between 0 and 20',
  },
  NumHAcceptors: {
    min: 0,
    max: 30,
    required: true,
    message: 'Number of H-acceptors must be between 0 and 30',
  },
  NumRotatableBonds: {
    min: 0,
    max: 50,
    required: true,
    message: 'Rotatable bonds must be between 0 and 50',
  },
  AromaticRings: {
    min: 0,
    max: 10,
    required: true,
    message: 'Aromatic rings must be between 0 and 10',
  },
  NumAromaticRings: {
    min: 0,
    max: 10,
    required: true,
    message: 'Aromatic rings must be between 0 and 10',
  },
  TPSA: {
    min: 0,
    max: 400,
    required: true,
    message: 'TPSA must be between 0 and 400',
  },
  NumHeteroatoms: {
    min: 0,
    max: 50,
    required: true,
    message: 'Heteroatoms must be between 0 and 50',
  },
  NumAromaticAtoms: {
    min: 0,
    max: 100,
    required: true,
    message: 'Aromatic atoms must be between 0 and 100',
  },
  NumSaturatedRings: {
    min: 0,
    max: 10,
    required: true,
    message: 'Saturated rings must be between 0 and 10',
  },
  NumAliphaticRings: {
    min: 0,
    max: 10,
    required: true,
    message: 'Aliphatic rings must be between 0 and 10',
  },
  RingCount: {
    min: 0,
    max: 10,
    required: true,
    message: 'Ring count must be between 0 and 10',
  },
  NumRings: {
    min: 0,
    max: 10,
    required: true,
    message: 'Number of rings must be between 0 and 10',
  },
  FractionCsp3: {
    min: 0,
    max: 1,
    required: true,
    message: 'Fraction Csp3 must be between 0 and 1',
  },
  FractionCSP3: {
    min: 0,
    max: 1,
    required: true,
    message: 'Fraction CSP3 must be between 0 and 1',
  },
  NumAromaticCarbocycles: {
    min: 0,
    max: 10,
    required: true,
    message: 'Aromatic carbocycles must be between 0 and 10',
  },
  NumSaturatedCarbocycles: {
    min: 0,
    max: 10,
    required: true,
    message: 'Saturated carbocycles must be between 0 and 10',
  },
  MolarRefractivity: {
    min: 0,
    max: 500,
    required: true,
    message: 'Molar refractivity must be between 0 and 500',
  },
  BertzCT: {
    min: 0,
    max: 5000,
    required: true,
    message: 'Bertz CT must be between 0 and 5000',
  },
  HeavyAtomCount: {
    min: 0,
    max: 150,
    required: true,
    message: 'Heavy atom count must be between 0 and 150',
  },
}

export const validateInput = (name: string, value: number): { isValid: boolean; message?: string } => {
  const rule = validationRules[name]
  if (!rule) return { isValid: true }

  if (rule.required && (value === null || value === undefined || isNaN(value))) {
    return { isValid: false, message: 'This field is required' }
  }

  if (rule.min !== undefined && value < rule.min) {
    return { isValid: false, message: rule.message || `Value must be at least ${rule.min}` }
  }

  if (rule.max !== undefined && value > rule.max) {
    return { isValid: false, message: rule.message || `Value must be at most ${rule.max}` }
  }

  return { isValid: true }
}

export const descriptorTooltips: Record<string, string> = {
  MolecularWeight: 'Mass of the molecule in Daltons (g/mol). Larger molecules tend to be less bioavailable.',
  LogP: 'Lipophilicity - measure of how well the compound dissolves in fat vs water. Higher values = more lipophilic (fat-soluble).',
  NumHDonors: 'Number of hydrogen bond donors (N-H, O-H groups). Important for biological interactions.',
  NumHAcceptors: 'Number of hydrogen bond acceptors (N, O atoms). Influences binding to proteins.',
  NumRotatableBonds: 'Number of bonds that can rotate freely. Higher values = more flexible molecule.',
  AromaticRings: 'Number of aromatic ring systems. Aromatic rings are stable, flat structures.',
  NumAromaticRings: 'Total count of aromatic rings in the molecule.',
  TPSA: 'Topological Polar Surface Area (Ų). Predicts membrane permeability and oral bioavailability.',
  NumHeteroatoms: 'Number of non-carbon/hydrogen atoms (N, O, S, P, etc.). Affects reactivity.',
  NumAromaticAtoms: 'Count of atoms in aromatic systems.',
  NumSaturatedRings: 'Number of rings with no double bonds.',
  NumAliphaticRings: 'Number of non-aromatic rings.',
  RingCount: 'Total number of ring systems in the molecule.',
  NumRings: 'Total count of all rings.',
  FractionCsp3: 'Fraction of carbon atoms with sp³ hybridization. Indicates 3D complexity.',
  FractionCSP3: 'Fraction of sp³-hybridized carbons. Higher = more 3D character.',
  NumAromaticCarbocycles: 'Number of carbocyclic aromatic rings (all-carbon aromatic rings).',
  NumSaturatedCarbocycles: 'Number of saturated carbocyclic rings.',
  MolarRefractivity: 'Measure of the volume occupied by atoms. Related to London dispersion forces.',
  BertzCT: 'Bertz complexity index. Measures overall molecular complexity.',
  HeavyAtomCount: 'Total number of non-hydrogen atoms. Basic measure of molecular size.',
}
