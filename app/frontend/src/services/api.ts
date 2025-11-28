import axios from 'axios'

// Use port 8001 (our fresh, working port!)
const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8001'

const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
})

export interface PredictionInput {
  source: string
  year: number
  toxicity_type: string
  herbicide: number
  fungicide: number
  insecticide: number
  other_agrochemical: number
  MolecularWeight: number
  LogP: number
  NumHDonors: number
  NumHAcceptors: number
  NumRotatableBonds: number
  AromaticRings: number
  NumAromaticRings: number
  TPSA: number
  NumHeteroatoms: number
  NumAromaticAtoms: number
  NumSaturatedRings: number
  NumAliphaticRings: number
  RingCount: number
  NumRings: number
  FractionCsp3: number
  FractionCSP3: number
  NumAromaticCarbocycles: number
  NumSaturatedCarbocycles: number
  MolarRefractivity: number
  BertzCT: number
  HeavyAtomCount: number
}

export interface PredictionResult {
  prediction: number
  label_text?: string
  prediction_label?: string
  probability_toxic?: number
  probability_non_toxic?: number
  probabilities?: {
    toxic: number
    non_toxic: number
  }
  confidence: number
  timestamp: string
}

export interface ModelInfo {
  model_name: string
  version: string
  algorithm: string
  metrics?: {
    test_accuracy?: number
    test_f1?: number
    test_roc_auc?: number
  }
}

export const predictToxicity = async (input: PredictionInput): Promise<PredictionResult> => {
  const response = await api.post<PredictionResult>('/predict', input)
  return response.data
}

export const getModelInfo = async (): Promise<ModelInfo> => {
  const response = await api.get<ModelInfo>('/model/info')
  return response.data
}

export const checkHealth = async (): Promise<{ status: string }> => {
  const response = await api.get('/health')
  return response.data
}

export default api
