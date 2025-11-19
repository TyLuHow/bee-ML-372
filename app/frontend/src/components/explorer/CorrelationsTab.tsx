import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { CorrelationsResponse } from '../../types/explorer'
import { HeatmapChart } from '../charts/HeatmapChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'

export const CorrelationsTab: React.FC = () => {
  const [data, setData] = useState<CorrelationsResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [selectedCell, setSelectedCell] = useState<{
    row: number
    col: number
    value: number
  } | null>(null)

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getCorrelations()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch correlations data')
    } finally {
      setLoading(false)
    }
  }

  const handleCellClick = (row: number, col: number, value: number) => {
    setSelectedCell({ row, col, value })
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

  // Sort top correlations by absolute value
  const sortedCorrelations = [...data.top_correlations].sort(
    (a, b) => Math.abs(b.correlation) - Math.abs(a.correlation)
  )

  return (
    <div className="space-y-6">
      {/* Introduction */}
      <Card>
        <CardContent className="p-6">
          <p className="text-text-secondary">
            Feature correlation analysis reveals relationships between molecular descriptors.
            Understanding these correlations helps identify redundant features and multicollinearity
            in predictive models.
          </p>
        </CardContent>
      </Card>

      {/* Correlation Heatmap */}
      <Card>
        <CardHeader>
          <CardTitle>Correlation Matrix Heatmap</CardTitle>
        </CardHeader>
        <CardContent>
          <HeatmapChart
            matrix={data.correlation_matrix}
            features={data.features}
            colorScale="diverging"
            onCellClick={handleCellClick}
          />
          {selectedCell && (
            <div className="mt-4 p-4 bg-honey-light/20 border border-bee-yellow/30 rounded-lg">
              <p className="font-semibold text-text-primary mb-2">Selected Correlation</p>
              <p className="text-sm text-text-secondary">
                <strong>{data.features[selectedCell.row]}</strong> vs{' '}
                <strong>{data.features[selectedCell.col]}</strong>:{' '}
                <span className={selectedCell.value > 0 ? 'text-accent-toxic' : 'text-accent-info'}>
                  {selectedCell.value.toFixed(3)}
                </span>
              </p>
            </div>
          )}
          <div className="mt-4 p-4 bg-bg-secondary rounded-lg">
            <p className="text-sm text-text-secondary">
              <strong>How to read:</strong> Red indicates positive correlation, blue indicates negative
              correlation, and white indicates no correlation. Click on any cell to see detailed values.
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Top Correlations List */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Positive Correlations */}
        <Card>
          <CardHeader>
            <CardTitle>Top Positive Correlations</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              {sortedCorrelations
                .filter(item => item.correlation > 0)
                .slice(0, 10)
                .map((item, index) => (
                  <div
                    key={index}
                    className="flex items-center justify-between p-3 bg-bg-secondary rounded-lg hover:bg-bg-tertiary transition-colors"
                  >
                    <div className="flex-1">
                      <p className="text-sm font-medium text-text-primary">
                        {item.feature1} ↔ {item.feature2}
                      </p>
                    </div>
                    <div className="flex items-center gap-3">
                      <span className="text-sm font-bold text-accent-toxic">
                        {item.correlation.toFixed(3)}
                      </span>
                      <span className="text-xs text-text-tertiary">
                        p={item.p_value < 0.001 ? '<0.001' : item.p_value.toFixed(3)}
                      </span>
                    </div>
                  </div>
                ))}
            </div>
          </CardContent>
        </Card>

        {/* Negative Correlations */}
        <Card>
          <CardHeader>
            <CardTitle>Top Negative Correlations</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              {sortedCorrelations
                .filter(item => item.correlation < 0)
                .slice(0, 10)
                .map((item, index) => (
                  <div
                    key={index}
                    className="flex items-center justify-between p-3 bg-bg-secondary rounded-lg hover:bg-bg-tertiary transition-colors"
                  >
                    <div className="flex-1">
                      <p className="text-sm font-medium text-text-primary">
                        {item.feature1} ↔ {item.feature2}
                      </p>
                    </div>
                    <div className="flex items-center gap-3">
                      <span className="text-sm font-bold text-accent-info">
                        {item.correlation.toFixed(3)}
                      </span>
                      <span className="text-xs text-text-tertiary">
                        p={item.p_value < 0.001 ? '<0.001' : item.p_value.toFixed(3)}
                      </span>
                    </div>
                  </div>
                ))}
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Full Correlation Table */}
      <Card>
        <CardHeader>
          <CardTitle>Top 20 Correlations (by Absolute Value)</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b-2 border-border-default">
                  <th className="text-left py-3 px-4 font-medium text-text-secondary">Rank</th>
                  <th className="text-left py-3 px-4 font-medium text-text-secondary">Feature 1</th>
                  <th className="text-left py-3 px-4 font-medium text-text-secondary">Feature 2</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Correlation</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">P-value</th>
                </tr>
              </thead>
              <tbody>
                {sortedCorrelations.slice(0, 20).map((item, index) => (
                  <tr key={index} className="border-b border-border-default hover:bg-bg-secondary">
                    <td className="py-3 px-4 text-text-tertiary">{index + 1}</td>
                    <td className="py-3 px-4 font-medium text-text-primary">{item.feature1}</td>
                    <td className="py-3 px-4 font-medium text-text-primary">{item.feature2}</td>
                    <td
                      className={`text-right py-3 px-4 font-bold ${
                        item.correlation > 0 ? 'text-accent-toxic' : 'text-accent-info'
                      }`}
                    >
                      {item.correlation.toFixed(3)}
                    </td>
                    <td className="text-right py-3 px-4 text-text-secondary">
                      {item.p_value < 0.001 ? '<0.001' : item.p_value.toFixed(3)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>

      {/* Interpretation Guide */}
      <Card>
        <CardHeader>
          <CardTitle>Interpretation Guide</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <h4 className="font-semibold text-text-primary mb-3">Correlation Strength</h4>
              <ul className="space-y-2 text-sm text-text-secondary">
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>|r| &gt; 0.7: Strong correlation</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>0.4 &lt; |r| &lt; 0.7: Moderate correlation</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>0.2 &lt; |r| &lt; 0.4: Weak correlation</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>|r| &lt; 0.2: Very weak or no correlation</span>
                </li>
              </ul>
            </div>
            <div>
              <h4 className="font-semibold text-text-primary mb-3">Practical Implications</h4>
              <ul className="space-y-2 text-sm text-text-secondary">
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>High correlations indicate redundant information</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>May require feature selection or dimensionality reduction</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Multicollinearity can affect model interpretation</span>
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-bee-yellow mt-1">•</span>
                  <span>Some algorithms handle correlations better than others</span>
                </li>
              </ul>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
