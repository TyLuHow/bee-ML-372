// API Response Types for Data Explorer

export interface OverviewResponse {
  total_compounds: number
  temporal_range: {
    min_year: number
    max_year: number
  }
  toxic_percentage: number
  data_sources: {
    name: string
    count: number
  }[]
  chemical_types: {
    type: string
    count: number
  }[]
  exposure_types: {
    type: string
    count: number
  }[]
}

export interface DistributionData {
  bin_center: number
  toxic: number
  non_toxic: number
}

export interface DescriptorStats {
  descriptor: string
  mean: number
  std: number
  min: number
  max: number
  toxic_mean: number
  non_toxic_mean: number
}

export interface MolecularDiversityResponse {
  descriptors: {
    [key: string]: {
      distribution: DistributionData[]
      stats: DescriptorStats
    }
  }
}

export interface ChemicalSpacePoint {
  compound_id: string
  compound_name: string
  x: number
  y: number
  toxicity: number
  year?: number
  chemical_type?: string
  logp?: number
  molecular_weight?: number
}

export interface ChemicalSpaceResponse {
  pca: {
    points: ChemicalSpacePoint[]
    variance_explained: number[]
  }
  tsne: {
    points: ChemicalSpacePoint[]
  }
}

export interface TemporalDataPoint {
  decade: string
  year: number
  count: number
  toxic_count: number
  toxicity_rate: number
}

export interface TemporalTrendsResponse {
  temporal_data: TemporalDataPoint[]
  rolling_average: {
    year: number
    toxicity_rate: number
  }[]
  mann_kendall: {
    tau: number
    p_value: number
    trend: string
  }
}

export interface ToxicophoreData {
  pattern: string
  pattern_name: string
  toxic_prevalence: number
  non_toxic_prevalence: number
  enrichment_ratio: number
  p_value: number
  ci_lower: number
  ci_upper: number
  toxic_count: number
  non_toxic_count: number
  toxicity_rate: number
}

export interface ToxicophoresResponse {
  toxicophores: ToxicophoreData[]
}

export interface CorrelationData {
  feature1: string
  feature2: string
  correlation: number
  p_value: number
}

export interface CorrelationsResponse {
  correlation_matrix: number[][]
  features: string[]
  top_correlations: CorrelationData[]
}

export interface PropertyScatterPoint {
  compound_id: string
  compound_name: string
  x: number
  y: number
  toxicity: number
}

export interface PropertyDistributionsResponse {
  scatter_plots: {
    logp_vs_mw: {
      data: PropertyScatterPoint[]
      x_label: string
      y_label: string
    }
    tpsa_vs_hdonors: {
      data: PropertyScatterPoint[]
      x_label: string
      y_label: string
    }
    rotatable_vs_csp3: {
      data: PropertyScatterPoint[]
      x_label: string
      y_label: string
    }
  }
}

// Chart data types
export interface PieChartData {
  name: string
  value: number
  color?: string
}

export interface BarChartData {
  name: string
  value: number
  error?: number
}

export interface LineChartData {
  x: number | string
  y: number
  label?: string
}

export interface ScatterChartData {
  x: number
  y: number
  label?: string
  color?: string
  [key: string]: any
}
