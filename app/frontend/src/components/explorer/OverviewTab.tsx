import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { OverviewResponse, PieChartData } from '../../types/explorer'
import { StatsCard } from '../charts/StatsCard'
import { PieChart } from '../charts/PieChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'

export const OverviewTab: React.FC = () => {
  const [data, setData] = useState<OverviewResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getOverview()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch overview data')
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="space-y-6">
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {Array.from({ length: 6 }).map((_, i) => (
            <Card key={i} className="animate-pulse">
              <CardContent className="p-6">
                <div className="h-4 bg-gray-200 rounded w-1/2 mb-4"></div>
                <div className="h-8 bg-gray-200 rounded w-3/4"></div>
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

  const chemicalTypesData: PieChartData[] = data.chemical_types.map(item => ({
    name: item.type,
    value: item.count,
  }))

  const sourcesData: PieChartData[] = data.data_sources.map(item => ({
    name: item.name,
    value: item.count,
  }))

  return (
    <div className="space-y-8">
      {/* Key Statistics */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
        <StatsCard
          title="Total Compounds"
          value={data.total_compounds.toLocaleString()}
          subtitle="Unique chemical structures"
          color="text-bee-yellow"
        />
        <StatsCard
          title="Temporal Range"
          value={`${data.temporal_range.min_year} - ${data.temporal_range.max_year}`}
          subtitle={`${data.temporal_range.max_year - data.temporal_range.min_year} years of data`}
          color="text-chart-primary"
        />
        <StatsCard
          title="Toxic Compounds"
          value={`${data.toxic_percentage.toFixed(1)}%`}
          subtitle={`${Math.round(data.total_compounds * data.toxic_percentage / 100)} compounds`}
          color="text-accent-toxic"
        />
        <StatsCard
          title="Data Sources"
          value={data.data_sources.length}
          subtitle="Scientific databases"
          color="text-chart-tertiary"
        />
        <StatsCard
          title="Chemical Types"
          value={data.chemical_types.length}
          subtitle="Unique categories"
          color="text-chart-secondary"
        />
        <StatsCard
          title="Exposure Types"
          value={data.exposure_types.length}
          subtitle="Test conditions"
          color="text-accent-info"
        />
      </div>

      {/* Distribution Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <Card>
          <CardHeader>
            <CardTitle>Chemical Type Distribution</CardTitle>
          </CardHeader>
          <CardContent>
            <PieChart
              data={chemicalTypesData}
              height={350}
            />
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Data Source Distribution</CardTitle>
          </CardHeader>
          <CardContent>
            <PieChart
              data={sourcesData}
              height={350}
            />
          </CardContent>
        </Card>
      </div>

      {/* Additional Info */}
      <Card>
        <CardHeader>
          <CardTitle>Dataset Information</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="prose max-w-none">
            <p className="text-text-secondary">
              The ApisTox dataset is a comprehensive collection of {data.total_compounds.toLocaleString()} unique chemical
              compounds with their molecular descriptors and toxicity labels for honey bees (Apis mellifera).
              The dataset spans {data.temporal_range.max_year - data.temporal_range.min_year} years of research,
              from {data.temporal_range.min_year} to {data.temporal_range.max_year}, and includes data from{' '}
              {data.data_sources.length} different scientific sources.
            </p>
            <p className="text-text-secondary mt-4">
              Approximately {data.toxic_percentage.toFixed(1)}% of the compounds in the dataset are classified as
              toxic to honey bees, making this a valuable resource for predictive modeling and understanding
              chemical-biological interactions in pollinator health research.
            </p>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
