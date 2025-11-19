import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { ChemicalSpaceResponse, ScatterChartData } from '../../types/explorer'
import { ScatterChart } from '../charts/ScatterChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'

type Method = 'pca' | 'tsne'
type ColorBy = 'toxicity' | 'year' | 'chemical_type' | 'logp'

export const ChemicalSpaceTab: React.FC = () => {
  const [data, setData] = useState<ChemicalSpaceResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [method, setMethod] = useState<Method>('pca')
  const [colorBy, setColorBy] = useState<ColorBy>('toxicity')

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getChemicalSpace()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch chemical space data')
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

  const points = method === 'pca' ? data.pca.points : data.tsne.points

  const chartData: ScatterChartData[] = points.map(point => {
    let color = undefined
    if (colorBy === 'toxicity') {
      color = point.toxicity === 1 ? 'hsl(0, 84%, 60%)' : 'hsl(142, 76%, 36%)'
    } else if (colorBy === 'logp' && point.logp !== undefined) {
      // Color by LogP quartiles
      const logp = point.logp
      if (logp < 2) color = 'hsl(217, 91%, 60%)'
      else if (logp < 4) color = 'hsl(188, 100%, 35%)'
      else if (logp < 6) color = 'hsl(31, 100%, 50%)'
      else color = 'hsl(0, 84%, 60%)'
    }

    return {
      x: point.x,
      y: point.y,
      toxicity: point.toxicity,
      label: point.compound_name,
      compound_id: point.compound_id,
      compound_name: point.compound_name,
      color,
    }
  })

  return (
    <div className="space-y-6">
      {/* Introduction */}
      <Card>
        <CardContent className="p-6">
          <p className="text-text-secondary">
            Chemical space visualization using dimensionality reduction techniques to map high-dimensional
            molecular descriptors onto 2D space. This reveals clustering patterns and structural similarities
            among toxic and non-toxic compounds.
          </p>
        </CardContent>
      </Card>

      {/* Controls */}
      <Card>
        <CardContent className="p-6">
          <div className="flex flex-wrap gap-4 items-center">
            {/* Method Toggle */}
            <div className="flex items-center gap-2">
              <span className="text-sm font-medium text-text-secondary">Method:</span>
              <div className="flex gap-2">
                <Button
                  variant={method === 'pca' ? 'primary' : 'outline'}
                  size="sm"
                  onClick={() => setMethod('pca')}
                >
                  PCA
                </Button>
                <Button
                  variant={method === 'tsne' ? 'primary' : 'outline'}
                  size="sm"
                  onClick={() => setMethod('tsne')}
                >
                  t-SNE
                </Button>
              </div>
            </div>

            {/* Color By Selector */}
            <div className="flex items-center gap-2">
              <span className="text-sm font-medium text-text-secondary">Color by:</span>
              <select
                value={colorBy}
                onChange={(e) => setColorBy(e.target.value as ColorBy)}
                className="px-3 py-1.5 border border-border-default rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-bee-yellow"
              >
                <option value="toxicity">Toxicity</option>
                <option value="year">Year</option>
                <option value="chemical_type">Chemical Type</option>
                <option value="logp">LogP Quartiles</option>
              </select>
            </div>

            {/* Variance Explained (PCA only) */}
            {method === 'pca' && data.pca.variance_explained && (
              <div className="ml-auto text-sm text-text-secondary">
                Variance explained: PC1: {(data.pca.variance_explained[0] * 100).toFixed(1)}%,
                PC2: {(data.pca.variance_explained[1] * 100).toFixed(1)}%
              </div>
            )}
          </div>
        </CardContent>
      </Card>

      {/* Scatter Plot */}
      <Card>
        <CardHeader>
          <CardTitle>
            {method === 'pca' ? 'Principal Component Analysis' : 't-SNE Projection'}
          </CardTitle>
        </CardHeader>
        <CardContent>
          <ScatterChart
            data={chartData}
            xLabel={method === 'pca' ? 'PC1' : 't-SNE 1'}
            yLabel={method === 'pca' ? 'PC2' : 't-SNE 2'}
            colorBy={colorBy === 'toxicity' ? 'toxicity' : 'category'}
            height={500}
            showLegend={colorBy === 'toxicity'}
          />
        </CardContent>
      </Card>

      {/* Method Information */}
      <Card>
        <CardHeader>
          <CardTitle>About {method === 'pca' ? 'PCA' : 't-SNE'}</CardTitle>
        </CardHeader>
        <CardContent>
          {method === 'pca' ? (
            <div className="text-text-secondary space-y-2">
              <p>
                <strong>Principal Component Analysis (PCA)</strong> is a linear dimensionality reduction
                technique that identifies the directions of maximum variance in the data. The first two
                principal components (PC1 and PC2) capture the most significant variations in molecular
                structure.
              </p>
              <p>
                The variance explained indicates how much of the original data's variability is preserved
                in the 2D projection. Higher percentages indicate better representation of the original
                high-dimensional space.
              </p>
            </div>
          ) : (
            <div className="text-text-secondary space-y-2">
              <p>
                <strong>t-Distributed Stochastic Neighbor Embedding (t-SNE)</strong> is a non-linear
                dimensionality reduction technique that excels at preserving local structure. It groups
                similar molecules together while separating dissimilar ones.
              </p>
              <p>
                t-SNE is particularly useful for identifying clusters and outliers in the chemical space.
                Note that distances between clusters may not be meaningful, unlike in PCA.
              </p>
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  )
}
