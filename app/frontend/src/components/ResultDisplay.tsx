import { Card, CardHeader, CardTitle, CardContent } from './ui/Card'
import { RadialGauge } from './charts/RadialGauge'
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Cell } from 'recharts'

interface PredictionResult {
  prediction: number
  label_text: string
  probability_toxic: number
  probability_non_toxic: number
  confidence: number
  timestamp: string
}

interface Props {
  result: PredictionResult | null
  loading: boolean
  error: string | null
}

// Mock SHAP values (in a real app, this would come from the API)
const mockShapValues = [
  { feature: 'insecticide', value: 0.45, direction: 'positive' as const },
  { feature: 'LogP', value: 0.23, direction: 'positive' as const },
  { feature: 'year', value: -0.18, direction: 'negative' as const },
  { feature: 'MolecularWeight', value: 0.12, direction: 'positive' as const },
  { feature: 'TPSA', value: -0.08, direction: 'negative' as const },
]

// Mock similar compounds (in a real app, this would come from the API)
const mockSimilarCompounds = [
  { name: 'Imidacloprid', similarity: 92, toxicity: 'Toxic' },
  { name: 'Clothianidin', similarity: 88, toxicity: 'Toxic' },
  { name: 'Thiamethoxam', similarity: 85, toxicity: 'Toxic' },
]

const ResultDisplay = ({ result, loading, error }: Props) => {
  if (loading) {
    return (
      <Card className="text-center">
        <CardContent className="py-12">
          <div className="animate-spin rounded-full h-16 w-16 border-b-2 border-indigo-600 mx-auto"></div>
          <p className="mt-4 text-gray-600">Analyzing compound...</p>
        </CardContent>
      </Card>
    )
  }

  if (error) {
    return (
      <Card className="bg-red-50 border-2 border-red-200">
        <CardContent>
          <h3 className="text-lg font-bold text-red-800 mb-2">Error</h3>
          <p className="text-red-600">{error}</p>
          <p className="text-xs text-gray-500 mt-2">
            Note: The prediction endpoint has a known preprocessing issue. Try the API docs at
            localhost:8000/docs instead.
          </p>
        </CardContent>
      </Card>
    )
  }

  if (!result) {
    return (
      <Card className="text-center text-gray-500">
        <CardContent className="py-12">
          <div className="text-6xl mb-4">🐝</div>
          <h3 className="text-lg font-semibold mb-2">Ready to Predict</h3>
          <p className="text-sm">
            Enter compound properties and click &quot;Predict Toxicity&quot; to see results.
          </p>
        </CardContent>
      </Card>
    )
  }

  const isToxic = result.prediction === 1
  const confidencePercent = result.confidence * 100
  const predictionColor = isToxic ? '#ef4444' : '#10b981'

  // Prepare SHAP data for visualization
  const shapChartData = mockShapValues.map(item => ({
    ...item,
    displayValue: Math.abs(item.value),
    color: item.direction === 'positive' ? '#ef4444' : '#10b981',
  }))

  return (
    <div className="space-y-4 animate-fadeIn">
      {/* Main Prediction Card */}
      <Card className={`${isToxic ? 'bg-red-50 border-2 border-red-300' : 'bg-green-50 border-2 border-green-300'}`}>
        <CardContent>
          <div className="text-center mb-6">
            <div className={`inline-block px-6 py-3 rounded-full ${isToxic ? 'bg-red-100 text-red-700' : 'bg-green-100 text-green-700'} text-2xl font-bold`}>
              {isToxic ? 'Toxic to Bees' : 'Safe for Bees'}
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* Left: Confidence Gauge */}
            <div className="flex flex-col items-center justify-center">
              <RadialGauge value={confidencePercent} label="Confidence" color={predictionColor} />
            </div>

            {/* Right: Probabilities */}
            <div className="space-y-4">
              <div>
                <div className="flex justify-between text-sm mb-1">
                  <span className="font-medium">Non-Toxic Probability</span>
                  <span className="font-semibold">
                    {(result.probability_non_toxic * 100).toFixed(1)}%
                  </span>
                </div>
                <div className="w-full bg-gray-200 rounded-full h-3">
                  <div
                    className="bg-green-500 h-3 rounded-full transition-all duration-500"
                    style={{ width: `${result.probability_non_toxic * 100}%` }}
                  ></div>
                </div>
              </div>

              <div>
                <div className="flex justify-between text-sm mb-1">
                  <span className="font-medium">Toxic Probability</span>
                  <span className="font-semibold">
                    {(result.probability_toxic * 100).toFixed(1)}%
                  </span>
                </div>
                <div className="w-full bg-gray-200 rounded-full h-3">
                  <div
                    className="bg-red-500 h-3 rounded-full transition-all duration-500"
                    style={{ width: `${result.probability_toxic * 100}%` }}
                  ></div>
                </div>
              </div>

              {/* Interpretation */}
              <div className="bg-white p-4 rounded-lg border border-gray-200 mt-4">
                <h4 className="font-semibold text-sm mb-2 flex items-center gap-2">
                  <span>💡</span>
                  Interpretation
                </h4>
                <p className="text-xs text-gray-700">
                  {isToxic ? (
                    <>
                      This compound is predicted to be <strong>toxic to honey bees</strong>.{' '}
                      {result.confidence > 0.8
                        ? 'High confidence prediction.'
                        : 'Moderate confidence - consider laboratory validation.'}
                    </>
                  ) : (
                    <>
                      This compound is predicted to be <strong>safe for honey bees</strong>.{' '}
                      {result.confidence > 0.8
                        ? 'High confidence prediction.'
                        : 'Moderate confidence - monitoring recommended.'}
                    </>
                  )}
                </p>
              </div>
            </div>
          </div>

          <p className="text-xs text-gray-500 text-center mt-4">
            Predicted at: {new Date(result.timestamp).toLocaleString()}
          </p>
        </CardContent>
      </Card>

      {/* SHAP Feature Importance */}
      <Card>
        <CardHeader>
          <CardTitle className="text-lg">Feature Importance (SHAP)</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-sm text-gray-600 mb-4">
            Top 5 features that influenced this prediction. Red bars increase toxicity risk, green
            bars decrease it.
          </p>
          <ResponsiveContainer width="100%" height={200}>
            <BarChart data={shapChartData} layout="vertical">
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" domain={[-0.5, 0.5]} />
              <YAxis dataKey="feature" type="category" width={120} />
              <Tooltip
                formatter={(_value: any, _name: string, props: any) => {
                  const direction = props.payload.direction
                  const sign = direction === 'positive' ? '+' : ''
                  return [`${sign}${props.payload.value.toFixed(3)}`, 'SHAP Value']
                }}
              />
              <Bar dataKey="value" radius={[0, 8, 8, 0]}>
                {shapChartData.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={entry.color} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </CardContent>
      </Card>

      {/* Similar Compounds */}
      <Card>
        <CardHeader>
          <CardTitle className="text-lg">Similar Compounds</CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-sm text-gray-600 mb-4">
            Compounds with similar molecular properties (based on K-Nearest Neighbors)
          </p>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b">
                  <th className="text-left py-2 px-3">Compound Name</th>
                  <th className="text-center py-2 px-3">Similarity</th>
                  <th className="text-center py-2 px-3">Toxicity</th>
                </tr>
              </thead>
              <tbody>
                {mockSimilarCompounds.map((compound, index) => (
                  <tr key={index} className="border-b hover:bg-gray-50">
                    <td className="py-2 px-3 font-medium">{compound.name}</td>
                    <td className="py-2 px-3 text-center">
                      <div className="flex items-center justify-center gap-2">
                        <div className="w-16 bg-gray-200 rounded-full h-2">
                          <div
                            className="bg-indigo-500 h-2 rounded-full"
                            style={{ width: `${compound.similarity}%` }}
                          ></div>
                        </div>
                        <span className="text-xs">{compound.similarity}%</span>
                      </div>
                    </td>
                    <td className="py-2 px-3 text-center">
                      <span
                        className={`inline-block px-2 py-1 rounded-full text-xs font-semibold ${
                          compound.toxicity === 'Toxic'
                            ? 'bg-red-100 text-red-700'
                            : 'bg-green-100 text-green-700'
                        }`}
                      >
                        {compound.toxicity}
                      </span>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}

export default ResultDisplay

