import {
  OverviewResponse,
  MolecularDiversityResponse,
  ChemicalSpaceResponse,
  TemporalTrendsResponse,
  ToxicophoresResponse,
  CorrelationsResponse,
  PropertyDistributionsResponse,
} from '../types/explorer'

const API_BASE_URL = import.meta.env.VITE_API_URL || '/api'

/**
 * Data Explorer API Service
 * Handles all data fetching for the Data Explorer page
 */
export const explorerApi = {
  /**
   * Fetch dataset overview statistics
   */
  getOverview: async (): Promise<OverviewResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/overview`)
    if (!response.ok) {
      throw new Error(`Failed to fetch overview: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Fetch molecular diversity distributions and statistics
   */
  getMolecularDiversity: async (): Promise<MolecularDiversityResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/molecular-diversity`)
    if (!response.ok) {
      throw new Error(`Failed to fetch molecular diversity: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Fetch chemical space coordinates (PCA and t-SNE)
   */
  getChemicalSpace: async (): Promise<ChemicalSpaceResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/chemical-space`)
    if (!response.ok) {
      throw new Error(`Failed to fetch chemical space: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Fetch temporal trends data
   */
  getTemporalTrends: async (): Promise<TemporalTrendsResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/temporal-trends`)
    if (!response.ok) {
      throw new Error(`Failed to fetch temporal trends: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Fetch toxicophore analysis data
   */
  getToxicophores: async (): Promise<ToxicophoresResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/toxicophores`)
    if (!response.ok) {
      throw new Error(`Failed to fetch toxicophores: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Fetch feature correlations
   */
  getCorrelations: async (): Promise<CorrelationsResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/correlations`)
    if (!response.ok) {
      throw new Error(`Failed to fetch correlations: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Fetch property distributions (scatter plots)
   */
  getPropertyDistributions: async (): Promise<PropertyDistributionsResponse> => {
    const response = await fetch(`${API_BASE_URL}/explorer/property-distributions`)
    if (!response.ok) {
      throw new Error(`Failed to fetch property distributions: ${response.statusText}`)
    }
    return response.json()
  },

  /**
   * Export data to CSV
   * @param endpoint - The explorer endpoint to export data from
   * @returns CSV blob for download
   */
  exportToCSV: async (endpoint: string): Promise<Blob> => {
    const response = await fetch(`${API_BASE_URL}/explorer/${endpoint}?format=csv`)
    if (!response.ok) {
      throw new Error(`Failed to export data: ${response.statusText}`)
    }
    return response.blob()
  },
}
