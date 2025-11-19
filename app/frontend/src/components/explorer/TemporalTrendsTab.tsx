import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { TemporalTrendsResponse, LineChartData, BarChartData } from '../../types/explorer'
import { LineChart } from '../charts/LineChart'
import { BarChart } from '../charts/BarChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'
import { Badge } from '../ui/Badge'

export const TemporalTrendsTab: React.FC = () => {
  const [data, setData] = useState<TemporalTrendsResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getTemporalTrends()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch temporal trends data')
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="space-y-6">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {Array.from({ length: 2 }).map((_, i) => (
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

  const toxicityRateData: LineChartData[] = data.temporal_data.map(point => ({
    x: point.decade,
    y: point.toxicity_rate * 100,
    label: 'Toxicity Rate',
  }))

  const countData: BarChartData[] = data.temporal_data.map(point => ({
    name: point.decade,
    value: point.count,
  }))

  const getTrendBadge = (trend: string) => {
    if (trend.toLowerCase().includes('increasing')) {
      return <Badge variant="toxic">Increasing</Badge>
    } else if (trend.toLowerCase().includes('decreasing')) {
      return <Badge variant="safe">Decreasing</Badge>
    } else {
      return <Badge variant="default">No Trend</Badge>
    }
  }

  const getSignificance = (pValue: number) => {
    if (pValue < 0.001) return '***'
    if (pValue < 0.01) return '**'
    if (pValue < 0.05) return '*'
    return 'ns'
  }

  return (
    <div className="space-y-6">
      {/* Introduction */}
      <Card>
        <CardContent className="p-6">
          <p className="text-text-secondary">
            Temporal analysis reveals how honey bee toxicity data has evolved over time,
            showing trends in research focus and chemical compound testing patterns across decades.
          </p>
        </CardContent>
      </Card>

      {/* Statistical Test Results */}
      <Card>
        <CardHeader>
          <CardTitle>Mann-Kendall Trend Test</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="p-4 bg-bg-secondary rounded-lg">
              <p className="text-sm text-text-tertiary mb-1">Trend Direction</p>
              <div className="flex items-center gap-2">
                {getTrendBadge(data.mann_kendall.trend)}
              </div>
            </div>
            <div className="p-4 bg-bg-secondary rounded-lg">
              <p className="text-sm text-text-tertiary mb-1">Kendall's Tau</p>
              <p className="text-2xl font-bold text-text-primary">
                {data.mann_kendall.tau.toFixed(3)}
              </p>
            </div>
            <div className="p-4 bg-bg-secondary rounded-lg">
              <p className="text-sm text-text-tertiary mb-1">P-value</p>
              <p className="text-2xl font-bold text-text-primary">
                {data.mann_kendall.p_value.toFixed(4)}
                <span className="text-sm ml-2">{getSignificance(data.mann_kendall.p_value)}</span>
              </p>
            </div>
          </div>
          <div className="mt-4 p-4 bg-honey-light/20 border border-bee-yellow/30 rounded-lg">
            <p className="text-sm text-text-secondary">
              <strong>Interpretation:</strong> {data.mann_kendall.trend}.
              {data.mann_kendall.p_value < 0.05
                ? ' The trend is statistically significant.'
                : ' The trend is not statistically significant.'}
            </p>
          </div>
        </CardContent>
      </Card>

      {/* Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Toxicity Rate Over Time */}
        <Card>
          <CardHeader>
            <CardTitle>Toxicity Rate by Decade</CardTitle>
          </CardHeader>
          <CardContent>
            <LineChart
              data={toxicityRateData}
              xLabel="Decade"
              yLabel="Toxicity Rate (%)"
              height={350}
              color="hsl(0, 84%, 60%)"
            />
          </CardContent>
        </Card>

        {/* Sample Count by Decade */}
        <Card>
          <CardHeader>
            <CardTitle>Sample Count by Decade</CardTitle>
          </CardHeader>
          <CardContent>
            <BarChart
              data={countData}
              xLabel="Decade"
              yLabel="Number of Compounds"
              height={350}
              color="hsl(262, 80%, 50%)"
            />
          </CardContent>
        </Card>
      </div>

      {/* Data Table */}
      <Card>
        <CardHeader>
          <CardTitle>Temporal Data Summary</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b-2 border-border-default">
                  <th className="text-left py-3 px-4 font-medium text-text-secondary">Decade</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Total Count</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Toxic Count</th>
                  <th className="text-right py-3 px-4 font-medium text-text-secondary">Toxicity Rate</th>
                </tr>
              </thead>
              <tbody>
                {data.temporal_data.map((row, index) => (
                  <tr key={index} className="border-b border-border-default hover:bg-bg-secondary">
                    <td className="py-3 px-4 font-medium text-text-primary">{row.decade}</td>
                    <td className="text-right py-3 px-4 text-text-primary">{row.count}</td>
                    <td className="text-right py-3 px-4 text-accent-toxic">{row.toxic_count}</td>
                    <td className="text-right py-3 px-4 text-text-primary">
                      {(row.toxicity_rate * 100).toFixed(1)}%
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>

      {/* Insights */}
      <Card>
        <CardHeader>
          <CardTitle>Key Insights</CardTitle>
        </CardHeader>
        <CardContent>
          <ul className="space-y-2 text-text-secondary">
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                The temporal analysis helps identify shifts in chemical testing priorities and
                regulatory focus over time.
              </span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                Mann-Kendall test detects monotonic trends while being robust to outliers and
                non-normal distributions.
              </span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                Variations in sample counts reflect changing research efforts and data availability
                across different time periods.
              </span>
            </li>
          </ul>
        </CardContent>
      </Card>
    </div>
  )
}
