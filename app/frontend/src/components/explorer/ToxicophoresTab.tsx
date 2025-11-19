import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { ToxicophoresResponse, BarChartData, ScatterChartData } from '../../types/explorer'
import { BarChart } from '../charts/BarChart'
import { ScatterChart } from '../charts/ScatterChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'

export const ToxicophoresTab: React.FC = () => {
  const [data, setData] = useState<ToxicophoresResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getToxicophores()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch toxicophores data')
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="space-y-6">
        <Card className="animate-pulse">
          <CardContent className="p-6">
            <div className="h-96 bg-gray-200 rounded"></div>
          </CardContent>
        </Card>
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

  const getSignificance = (pValue: number): string => {
    if (pValue < 0.001) return '***'
    if (pValue < 0.01) return '**'
    if (pValue < 0.05) return '*'
    return ''
  }

  // Sort by enrichment ratio
  const sortedToxicophores = [...data.toxicophores].sort(
    (a, b) => b.enrichment_ratio - a.enrichment_ratio
  )

  const enrichmentData: BarChartData[] = sortedToxicophores.slice(0, 15).map(item => ({
    name: item.pattern_name,
    value: item.enrichment_ratio,
    error: item.ci_upper - item.enrichment_ratio,
  }))

  const scatterData: ScatterChartData[] = data.toxicophores.map(item => {
    const totalCount = item.toxic_count + item.non_toxic_count
    return {
      x: item.toxic_prevalence * 100,
      y: item.toxicity_rate * 100,
      label: item.pattern_name,
      size: totalCount,
      toxicity: item.toxicity_rate > 0.5 ? 1 : 0,
    }
  })

  return (
    <div className="space-y-6">
      {/* Introduction */}
      <Card>
        <CardContent className="p-6">
          <p className="text-text-secondary">
            Toxicophores are structural alerts - molecular substructures associated with toxic activity.
            This analysis identifies chemical patterns that are significantly enriched in toxic compounds
            compared to non-toxic ones.
          </p>
        </CardContent>
      </Card>

      {/* Enrichment Bar Chart */}
      <Card>
        <CardHeader>
          <CardTitle>Top Toxicophores by Enrichment Ratio</CardTitle>
        </CardHeader>
        <CardContent>
          <BarChart
            data={enrichmentData}
            orientation="horizontal"
            color="hsl(0, 84%, 60%)"
            showError={true}
            xLabel="Enrichment Ratio"
            height={500}
          />
          <div className="mt-4 p-4 bg-bg-secondary rounded-lg">
            <p className="text-sm text-text-secondary">
              <strong>Enrichment Ratio:</strong> The fold-change in prevalence of a structural pattern
              in toxic vs non-toxic compounds. Values &gt; 1 indicate enrichment in toxic compounds.
              Error bars show 95% confidence intervals.
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Scatter Plot */}
      <Card>
        <CardHeader>
          <CardTitle>Prevalence vs Toxicity Rate</CardTitle>
        </CardHeader>
        <CardContent>
          <ScatterChart
            data={scatterData}
            xLabel="Prevalence in Toxic Compounds (%)"
            yLabel="Toxicity Rate (%)"
            colorBy="toxicity"
            height={450}
          />
          <div className="mt-4 p-4 bg-bg-secondary rounded-lg">
            <p className="text-sm text-text-secondary">
              <strong>Interpretation:</strong> Points in the upper right represent highly prevalent
              and highly predictive toxicophores. These are the most useful structural alerts for
              toxicity prediction.
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Top Toxicophores Table */}
      <Card>
        <CardHeader>
          <CardTitle>Top 10 Toxicophores</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b-2 border-border-default">
                  <th className="text-left py-3 px-4 font-medium text-text-secondary">Pattern</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Enrichment</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Toxic Prev.</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Non-Toxic Prev.</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">P-value</th>
                  <th className="text-center py-3 px-4 font-medium text-text-secondary">Sig.</th>
                </tr>
              </thead>
              <tbody>
                {sortedToxicophores.slice(0, 10).map((item, index) => (
                  <tr key={index} className="border-b border-border-default hover:bg-bg-secondary">
                    <td className="py-3 px-4 font-medium text-text-primary">{item.pattern_name}</td>
                    <td className="text-right py-3 px-4 text-accent-toxic font-semibold">
                      {item.enrichment_ratio.toFixed(2)}
                    </td>
                    <td className="text-right py-3 px-4 text-text-primary">
                      {(item.toxic_prevalence * 100).toFixed(1)}%
                    </td>
                    <td className="text-right py-3 px-4 text-text-primary">
                      {(item.non_toxic_prevalence * 100).toFixed(1)}%
                    </td>
                    <td className="text-right py-3 px-4 text-text-primary">
                      {item.p_value < 0.001 ? '<0.001' : item.p_value.toFixed(3)}
                    </td>
                    <td className="text-center py-3 px-4 text-bee-yellow font-bold">
                      {getSignificance(item.p_value)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="mt-4 text-xs text-text-tertiary">
            <p>Significance levels: *** p &lt; 0.001, ** p &lt; 0.01, * p &lt; 0.05</p>
          </div>
        </CardContent>
      </Card>

      {/* Statistical Details */}
      <Card>
        <CardHeader>
          <CardTitle>Statistical Details</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <h4 className="font-semibold text-text-primary mb-3">Methodology</h4>
              <ul className="space-y-2 text-sm text-text-secondary">
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Structural patterns identified using SMARTS substructure matching</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Enrichment ratios calculated as (toxic_prev / non_toxic_prev)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Statistical significance determined by Fisher's exact test</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>95% confidence intervals computed using bootstrap resampling</span>
                </li>
              </ul>
            </div>
            <div>
              <h4 className="font-semibold text-text-primary mb-3">Applications</h4>
              <ul className="space-y-2 text-sm text-text-secondary">
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Early screening for potentially toxic compounds</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Guiding safer chemical design (avoid toxicophores)</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Feature engineering for machine learning models</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Understanding mechanisms of action</span>
                </li>
              </ul>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
