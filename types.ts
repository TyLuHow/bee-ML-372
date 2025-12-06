export interface ChemicalData {
  name: string;
  smiles?: string;
  mw: number;
  logP: number;
  exposure: string;
  category: string;
}

export interface PredictionResult {
  toxicity: 'Toxic' | 'Safe' | 'Uncertain';
  confidence: number;
  explanation: string;
  recommendation: string;
}

export type TabView = 'engine' | 'scenarios' | 'science';

export interface Scenario {
  id: string;
  title: string;
  role: string;
  description: string;
  icon: string;
}