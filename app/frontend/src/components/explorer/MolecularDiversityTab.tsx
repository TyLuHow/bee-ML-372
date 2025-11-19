import React, { useState, useEffect } from 'react'
import { explorerApi } from '../../services/explorerApi'
import { MolecularDiversityResponse } from '../../types/explorer'
import { DistributionChart } from '../charts/DistributionChart'
import { Card, CardHeader, CardTitle, CardContent } from '../ui/Card'
import { Alert } from '../ui/Alert'
import { Button } from '../ui/Button'

export const MolecularDiversityTab: React.FC = () => {
  const [data, setData] = useState<MolecularDiversityResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    fetchData()
  }, [])

  const fetchData = async () => {
    try {
      setLoading(true)
      setError(null)
      const result = await explorerApi.getMolecularDiversity()
      setData(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to fetch molecular diversity data')
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="space-y-6">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {Array.from({ length: 6 }).map((_, i) => (
            <Card key={i} className="animate-pulse">
              <CardContent className="p-6">
                <div className="h-64 bg-gray-200 rounded"></div>
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

  const descriptorNames = Object.keys(data.descriptors)

  return (
    <div className="space-y-8">
      {/* Introduction */}
      <Card>
        <CardContent className="p-6">
          <p className="text-text-secondary">
            These distributions show the molecular diversity across key chemical descriptors.
            Each histogram compares toxic (red) vs non-toxic (green) compounds to identify
            structural features associated with honey bee toxicity.
          </p>
        </CardContent>
      </Card>

      {/* Distribution Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {descriptorNames.map(descriptorName => {
          const descriptor = data.descriptors[descriptorName]
          return (
            <Card key={descriptorName}>
              <CardHeader>
                <CardTitle>{descriptorName}</CardTitle>
              </CardHeader>
              <CardContent>
                <DistributionChart
                  data={descriptor.distribution}
                  xLabel={descriptorName}
                  yLabel="Count"
                  height={300}
                />
                {/* Statistics Table */}
                <div className="mt-4 overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b border-border-default">
                        <th className="text-left py-2 text-text-secondary font-medium">Statistic</th>
                        <th className="text-right py-2 text-text-secondary font-medium">All</th>
                        <th className="text-right py-2 text-accent-toxic font-medium">Toxic</th>
                        <th className="text-right py-2 text-accent-safe font-medium">Non-Toxic</th>
                      </tr>
                    </thead>
                    <tbody className="text-text-primary">
                      <tr className="border-b border-border-default">
                        <td className="py-2">Mean</td>
                        <td className="text-right">{descriptor.stats.mean.toFixed(2)}</td>
                        <td className="text-right">{descriptor.stats.toxic_mean.toFixed(2)}</td>
                        <td className="text-right">{descriptor.stats.non_toxic_mean.toFixed(2)}</td>
                      </tr>
                      <tr className="border-b border-border-default">
                        <td className="py-2">Std Dev</td>
                        <td className="text-right">{descriptor.stats.std.toFixed(2)}</td>
                        <td className="text-right" colSpan={2}>-</td>
                      </tr>
                      <tr>
                        <td className="py-2">Range</td>
                        <td className="text-right" colSpan={3}>
                          {descriptor.stats.min.toFixed(2)} - {descriptor.stats.max.toFixed(2)}
                        </td>
                      </tr>
                    </tbody>
                  </table>
                </div>
              </CardContent>
            </Card>
          )
        })}
      </div>

      {/* Key Insights */}
      <Card>
        <CardHeader>
          <CardTitle>Key Insights</CardTitle>
        </CardHeader>
        <CardContent>
          <ul className="space-y-2 text-text-secondary">
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                Molecular weight distributions reveal structural complexity differences between toxic and non-toxic compounds.
              </span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                LogP (lipophilicity) values indicate cell membrane permeability, a key factor in bioavailability.
              </span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                Polar surface area (TPSA) affects absorption and blood-brain barrier penetration.
              </span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-bee-yellow mt-1">•</span>
              <span>
                Hydrogen bond donors/acceptors influence molecular interactions with biological targets.
              </span>
            </li>
          </ul>
        </CardContent>
      </Card>
    </div>
  )
}
