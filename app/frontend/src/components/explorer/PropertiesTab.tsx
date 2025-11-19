import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { PropertyDistributionsResponse, ScatterChartData } from '../../types/explorer'
import { ScatterChart } from '../charts/ScatterChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'

export const PropertiesTab: React.FC = () => {
  const [data, setData] = useState<PropertyDistributionsResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getPropertyDistributions()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch property distributions data')
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="space-y-6">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {Array.from({ length: 3 }).map((_, i) => (
            <Card key={i} className="animate-pulse">
              <CardContent className="p-6">
                <div className="h-80 bg-gray-200 rounded"></div>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>
    )
  }

  if (error || !data) {
    return (
      <div className="space-y-4">
        <Alert variant="error">
          <p className="font-medium">Failed to load data</p>
          <p className="text-sm mt-1">{error || 'Unknown error occurred'}</p>
        </Alert>
        <Button onClick={fetchData} variant="primary">
          Retry
        </Button>
      </div>
    )
  }

  const logpMwData: ScatterChartData[] = data.scatter_plots.logp_vs_mw.data.map(point => ({
    x: point.x,
    y: point.y,
    toxicity: point.toxicity,
    compound_id: point.compound_id,
    compound_name: point.compound_name,
  }))

  const tpsaHdonorsData: ScatterChartData[] = data.scatter_plots.tpsa_vs_hdonors.data.map(point => ({
    x: point.x,
    y: point.y,
    toxicity: point.toxicity,
    compound_id: point.compound_id,
    compound_name: point.compound_name,
  }))

  const rotatableCsp3Data: ScatterChartData[] = data.scatter_plots.rotatable_vs_csp3.data.map(point => ({
    x: point.x,
    y: point.y,
    toxicity: point.toxicity,
    compound_id: point.compound_id,
    compound_name: point.compound_name,
  }))

  return (
    <div className="space-y-6">
      {/* Introduction */}
      <Card>
        <CardContent className="p-6">
          <p className="text-text-secondary">
            Pairwise property distributions reveal relationships between key molecular descriptors
            and their association with toxicity. These scatter plots help identify chemical space
            regions enriched with toxic or non-toxic compounds.
          </p>
        </CardContent>
      </Card>

      {/* Scatter Plots Grid */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* LogP vs Molecular Weight */}
        <Card>
          <CardHeader>
            <CardTitle>Lipophilicity vs Size</CardTitle>
          </CardHeader>
          <CardContent>
            <ScatterChart
              data={logpMwData}
              xLabel={data.scatter_plots.logp_vs_mw.x_label}
              yLabel={data.scatter_plots.logp_vs_mw.y_label}
              colorBy="toxicity"
              height={350}
            />
          </CardContent>
        </Card>

        {/* TPSA vs H-Donors */}
        <Card>
          <CardHeader>
            <CardTitle>Polarity vs H-Bond Donors</CardTitle>
          </CardHeader>
          <CardContent>
            <ScatterChart
              data={tpsaHdonorsData}
              xLabel={data.scatter_plots.tpsa_vs_hdonors.x_label}
              yLabel={data.scatter_plots.tpsa_vs_hdonors.y_label}
              colorBy="toxicity"
              height={350}
            />
          </CardContent>
        </Card>

        {/* Rotatable Bonds vs Fraction CSP3 */}
        <Card>
          <CardHeader>
            <CardTitle>Flexibility vs Saturation</CardTitle>
          </CardHeader>
          <CardContent>
            <ScatterChart
              data={rotatableCsp3Data}
              xLabel={data.scatter_plots.rotatable_vs_csp3.x_label}
              yLabel={data.scatter_plots.rotatable_vs_csp3.y_label}
              colorBy="toxicity"
              height={350}
            />
          </CardContent>
        </Card>
      </div>

      {/* Property Descriptions */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <Card>
          <CardHeader>
            <CardTitle>LogP & Molecular Weight</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3 text-sm text-text-secondary">
              <p>
                <strong className="text-text-primary">LogP (Lipophilicity):</strong> Measures how
                hydrophobic a molecule is. Higher values indicate greater fat solubility, affecting
                membrane permeability and bioaccumulation.
              </p>
              <p>
                <strong className="text-text-primary">Molecular Weight:</strong> The mass of the
                molecule. Larger molecules may have different absorption and distribution profiles.
              </p>
              <p className="text-xs text-text-tertiary mt-2">
                Lipinski's Rule: LogP &lt; 5, MW &lt; 500 for drug-like compounds
              </p>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>TPSA & H-Bond Donors</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3 text-sm text-text-secondary">
              <p>
                <strong className="text-text-primary">TPSA (Topological Polar Surface Area):</strong>{' '}
                Sum of surface areas of polar atoms. Affects cell membrane permeability and
                blood-brain barrier penetration.
              </p>
              <p>
                <strong className="text-text-primary">H-Bond Donors:</strong> Number of hydrogen
                atoms attached to electronegative atoms. Influences solubility and binding affinity.
              </p>
              <p className="text-xs text-text-tertiary mt-2">
                Optimal range: TPSA 20-130 Ų, H-donors ≤ 5
              </p>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Flexibility & Saturation</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3 text-sm text-text-secondary">
              <p>
                <strong className="text-text-primary">Rotatable Bonds:</strong> Number of bonds
                allowing free rotation. Higher values indicate greater molecular flexibility,
                affecting binding specificity.
              </p>
              <p>
                <strong className="text-text-primary">Fraction CSP3:</strong> Proportion of sp³
                hybridized carbons. Higher values indicate more saturated, 3D structure vs flat
                aromatic rings.
              </p>
              <p className="text-xs text-text-tertiary mt-2">
                Flexible molecules may have multiple conformations
              </p>
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Key Insights */}
      <Card>
        <CardHeader>
          <CardTitle>Key Insights from Property Distributions</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <h4 className="font-semibold text-text-primary mb-3">Toxicity Patterns</h4>
              <ul className="space-y-2 text-sm text-text-secondary">
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>
                    Toxic compounds may cluster in specific regions of chemical property space
                  </span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>
                    Extreme values (very high LogP, very large MW) may correlate with toxicity
                  </span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>
                    Violations of drug-like rules may indicate problematic compounds
                  </span>
                </li>
              </ul>
            </div>
            <div>
              <h4 className="font-semibold text-text-primary mb-3">Model Implications</h4>
              <ul className="space-y-2 text-sm text-text-secondary">
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>
                    These properties serve as important features in machine learning models
                  </span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>
                    Non-linear relationships suggest need for complex models (e.g., neural networks)
                  </span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>
                    Property combinations may be more informative than individual properties
                  </span>
                </li>
              </ul>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Export Option */}
      <Card>
        <CardContent className="p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="font-medium text-text-primary">Export Property Data</p>
              <p className="text-sm text-text-secondary mt-1">
                Download scatter plot data for further analysis in external tools
              </p>
            </div>
            <Button variant="outline" size="sm">
              Export CSV
            </Button>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
