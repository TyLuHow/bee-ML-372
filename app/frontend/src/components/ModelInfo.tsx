import { useEffect, useState } from 'react'
import { getModelInfo, ModelInfo as ModelInfoType } from '../services/api'
import { Card, CardHeader, CardTitle, CardContent } from './ui/Card'
import { Tabs, TabsList, TabsTrigger, TabsContent } from './ui/Tabs'
import { ConfusionMatrix } from './charts/ConfusionMatrix'
import {
  LineChart,
  Line,
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  Cell,
} from 'recharts'

// Mock data for confusion matrix
const confusionMatrixData = {
  truePositive: 89,
  falsePositive: 18,
  trueNegative: 84,
  falseNegative: 16,
}

// Mock ROC curve data
const rocCurveData = [
  { fpr: 0, tpr: 0 },
  { fpr: 0.05, tpr: 0.45 },
  { fpr: 0.1, tpr: 0.68 },
  { fpr: 0.15, tpr: 0.78 },
  { fpr: 0.2, tpr: 0.84 },
  { fpr: 0.25, tpr: 0.88 },
  { fpr: 0.3, tpr: 0.91 },
  { fpr: 0.4, tpr: 0.94 },
  { fpr: 0.5, tpr: 0.96 },
  { fpr: 0.7, tpr: 0.98 },
  { fpr: 1, tpr: 1 },
]

// Mock Precision-Recall curve data
const prCurveData = [
  { recall: 0, precision: 1 },
  { recall: 0.1, precision: 0.98 },
  { recall: 0.2, precision: 0.95 },
  { recall: 0.3, precision: 0.92 },
  { recall: 0.4, precision: 0.89 },
  { recall: 0.5, precision: 0.86 },
  { recall: 0.6, precision: 0.82 },
  { recall: 0.7, precision: 0.78 },
  { recall: 0.8, precision: 0.72 },
  { recall: 0.9, precision: 0.65 },
  { recall: 1, precision: 0.55 },
]

// Mock feature importance data
const featureImportanceData = [
  { feature: 'insecticide', importance: 0.285 },
  { feature: 'herbicide', importance: 0.182 },
  { feature: 'fungicide', importance: 0.158 },
  { feature: 'year', importance: 0.095 },
  { feature: 'LogP', importance: 0.072 },
  { feature: 'MolecularWeight', importance: 0.048 },
  { feature: 'TPSA', importance: 0.039 },
  { feature: 'NumHDonors', importance: 0.028 },
  { feature: 'NumHAcceptors', importance: 0.025 },
  { feature: 'NumRotatableBonds', importance: 0.019 },
  { feature: 'NumRings', importance: 0.016 },
  { feature: 'FractionCSP3', importance: 0.012 },
  { feature: 'MolarRefractivity', importance: 0.010 },
  { feature: 'BertzCT', importance: 0.007 },
  { feature: 'HeavyAtomCount', importance: 0.004 },
]

// Mock model comparison data
const modelComparisonData = [
  {
    model: 'XGBoost',
    accuracy: 0.836,
    precision: 0.832,
    recall: 0.848,
    f1: 0.840,
    rocAuc: 0.858,
    trainingTime: '2.3s',
  },
  {
    model: 'Random Forest',
    accuracy: 0.821,
    precision: 0.815,
    recall: 0.829,
    f1: 0.822,
    rocAuc: 0.843,
    trainingTime: '4.1s',
  },
  {
    model: 'Logistic Regression',
    accuracy: 0.789,
    precision: 0.782,
    recall: 0.798,
    f1: 0.790,
    rocAuc: 0.812,
    trainingTime: '0.5s',
  },
  {
    model: 'SVM',
    accuracy: 0.803,
    precision: 0.798,
    recall: 0.811,
    f1: 0.804,
    rocAuc: 0.825,
    trainingTime: '3.7s',
  },
  {
    model: 'Neural Network',
    accuracy: 0.812,
    precision: 0.806,
    recall: 0.821,
    f1: 0.813,
    rocAuc: 0.834,
    trainingTime: '8.2s',
  },
  {
    model: 'Naive Bayes',
    accuracy: 0.765,
    precision: 0.758,
    recall: 0.775,
    f1: 0.766,
    rocAuc: 0.789,
    trainingTime: '0.3s',
  },
]

const ModelInfo = () => {
  const [info, setInfo] = useState<ModelInfoType | null>(null)

  useEffect(() => {
    const fetchInfo = async () => {
      try {
        const data = await getModelInfo()
        setInfo(data)
      } catch (err) {
        // Silently fail - will use fallback data
        console.error('Failed to fetch model info:', err)
      }
    }
    fetchInfo()
  }, [])

  return (
    <Card>
      <CardHeader>
        <CardTitle>Model Information</CardTitle>
      </CardHeader>
      <CardContent>
        <Tabs defaultValue="overview">
          <TabsList>
            <TabsTrigger value="overview">Overview</TabsTrigger>
            <TabsTrigger value="performance">Performance</TabsTrigger>
            <TabsTrigger value="features">Features</TabsTrigger>
            <TabsTrigger value="comparison">Model Comparison</TabsTrigger>
          </TabsList>

          {/* Overview Tab */}
          <TabsContent value="overview">
            <div className="space-y-4">
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div className="bg-gradient-to-r from-indigo-50 to-purple-50 p-4 rounded-lg">
                  <h4 className="text-sm font-semibold text-gray-700 mb-2">Model Details</h4>
                  <div className="space-y-2 text-sm">
                    <div>
                      <span className="font-medium">Algorithm:</span>{' '}
                      {info?.algorithm || 'XGBoost'}
                    </div>
                    <div>
                      <span className="font-medium">Version:</span> {info?.version || '1.0.0'}
                    </div>
                    <div>
                      <span className="font-medium">Training Date:</span> January 2025
                    </div>
                  </div>
                </div>

                <div className="bg-gradient-to-r from-green-50 to-blue-50 p-4 rounded-lg">
                  <h4 className="text-sm font-semibold text-gray-700 mb-2">Dataset</h4>
                  <div className="space-y-2 text-sm">
                    <div>
                      <span className="font-medium">Total Compounds:</span> 1,035
                    </div>
                    <div>
                      <span className="font-medium">Training Set:</span> 828 (80%)
                    </div>
                    <div>
                      <span className="font-medium">Test Set:</span> 207 (20%)
                    </div>
                  </div>
                </div>

                <div className="bg-gradient-to-r from-yellow-50 to-orange-50 p-4 rounded-lg">
                  <h4 className="text-sm font-semibold text-gray-700 mb-2">Features</h4>
                  <div className="space-y-2 text-sm">
                    <div>
                      <span className="font-medium">Total Features:</span> 24
                    </div>
                    <div>
                      <span className="font-medium">Molecular Descriptors:</span> 15
                    </div>
                    <div>
                      <span className="font-medium">Chemical Flags:</span> 4
                    </div>
                  </div>
                </div>

                <div className="bg-gradient-to-r from-pink-50 to-red-50 p-4 rounded-lg">
                  <h4 className="text-sm font-semibold text-gray-700 mb-2">Performance</h4>
                  <div className="space-y-2 text-sm">
                    <div>
                      <span className="font-medium">Accuracy:</span>{' '}
                      {info?.metrics?.test_accuracy
                        ? (info.metrics.test_accuracy * 100).toFixed(1)
                        : '83.6'}
                      %
                    </div>
                    <div>
                      <span className="font-medium">ROC-AUC:</span>{' '}
                      {info?.metrics?.test_roc_auc
                        ? (info.metrics.test_roc_auc * 100).toFixed(1)
                        : '85.8'}
                      %
                    </div>
                    <div>
                      <span className="font-medium">F1 Score:</span>{' '}
                      {info?.metrics?.test_f1 ? (info.metrics.test_f1 * 100).toFixed(1) : '84.0'}%
                    </div>
                  </div>
                </div>
              </div>

              <div className="bg-gray-50 p-4 rounded-lg">
                <h4 className="font-semibold text-sm text-gray-700 mb-2">Top 5 Predictors</h4>
                <ol className="text-sm text-gray-600 space-y-1 list-decimal list-inside">
                  <li>Insecticide flag (28.5% importance)</li>
                  <li>Herbicide flag (18.2% importance)</li>
                  <li>Fungicide flag (15.8% importance)</li>
                  <li>Publication year (9.5% importance)</li>
                  <li>LogP - Lipophilicity (7.2% importance)</li>
                </ol>
              </div>

              <div className="text-center">
                <a
                  href="http://localhost:8000/docs"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="inline-block text-sm text-indigo-600 hover:text-indigo-800 underline"
                >
                  View API Documentation
                </a>
              </div>
            </div>
          </TabsContent>

          {/* Performance Tab */}
          <TabsContent value="performance">
            <div className="space-y-6">
              {/* Confusion Matrix */}
              <ConfusionMatrix {...confusionMatrixData} />

              {/* ROC and PR Curves */}
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <div>
                  <h4 className="text-lg font-semibold text-gray-800 mb-4 text-center">
                    ROC Curve (AUC = 0.858)
                  </h4>
                  <ResponsiveContainer width="100%" height={300}>
                    <LineChart data={rocCurveData}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis
                        dataKey="fpr"
                        label={{ value: 'False Positive Rate', position: 'insideBottom', offset: -5 }}
                      />
                      <YAxis label={{ value: 'True Positive Rate', angle: -90, position: 'insideLeft' }} />
                      <Tooltip />
                      <Legend />
                      <Line
                        type="monotone"
                        dataKey="tpr"
                        stroke="#8b5cf6"
                        strokeWidth={3}
                        name="ROC Curve"
                        dot={false}
                      />
                      <Line
                        type="monotone"
                        data={[
                          { fpr: 0, tpr: 0 },
                          { fpr: 1, tpr: 1 },
                        ]}
                        dataKey="tpr"
                        stroke="#9ca3af"
                        strokeDasharray="5 5"
                        strokeWidth={2}
                        name="Random Classifier"
                        dot={false}
                      />
                    </LineChart>
                  </ResponsiveContainer>
                </div>

                <div>
                  <h4 className="text-lg font-semibold text-gray-800 mb-4 text-center">
                    Precision-Recall Curve (AP = 0.864)
                  </h4>
                  <ResponsiveContainer width="100%" height={300}>
                    <LineChart data={prCurveData}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis
                        dataKey="recall"
                        label={{ value: 'Recall', position: 'insideBottom', offset: -5 }}
                      />
                      <YAxis label={{ value: 'Precision', angle: -90, position: 'insideLeft' }} />
                      <Tooltip />
                      <Legend />
                      <Line
                        type="monotone"
                        dataKey="precision"
                        stroke="#10b981"
                        strokeWidth={3}
                        name="PR Curve"
                        dot={false}
                      />
                    </LineChart>
                  </ResponsiveContainer>
                </div>
              </div>
            </div>
          </TabsContent>

          {/* Features Tab */}
          <TabsContent value="features">
            <div className="space-y-4">
              <div className="bg-gray-50 p-4 rounded-lg">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  Feature Importance (Top 15)
                </h4>
                <p className="text-sm text-gray-600 mb-4">
                  Features ranked by their contribution to model predictions
                </p>
              </div>

              <ResponsiveContainer width="100%" height={500}>
                <BarChart data={featureImportanceData} layout="vertical">
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis type="number" label={{ value: 'Importance', position: 'insideBottom', offset: -5 }} />
                  <YAxis dataKey="feature" type="category" width={150} />
                  <Tooltip formatter={(value: any) => `${(value * 100).toFixed(1)}%`} />
                  <Bar dataKey="importance" radius={[0, 8, 8, 0]}>
                    {featureImportanceData.map((_entry, index) => (
                      <Cell
                        key={`cell-${index}`}
                        fill={
                          index < 3
                            ? '#8b5cf6'
                            : index < 7
                            ? '#6366f1'
                            : index < 10
                            ? '#3b82f6'
                            : '#60a5fa'
                        }
                      />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>

              <div className="bg-yellow-50 p-4 rounded-lg border border-yellow-200">
                <h5 className="font-semibold text-sm text-gray-800 mb-2">Key Insights</h5>
                <ul className="text-sm text-gray-700 space-y-1 list-disc list-inside">
                  <li>Chemical type flags (insecticide, herbicide, fungicide) dominate predictions</li>
                  <li>Temporal factor (year) is significant - newer compounds tend to be more toxic</li>
                  <li>Lipophilicity (LogP) is the most important molecular descriptor</li>
                  <li>Molecular weight and polarity (TPSA) have moderate influence</li>
                </ul>
              </div>
            </div>
          </TabsContent>

          {/* Model Comparison Tab */}
          <TabsContent value="comparison">
            <div className="space-y-4">
              <div className="bg-gray-50 p-4 rounded-lg">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  Model Performance Comparison
                </h4>
                <p className="text-sm text-gray-600">
                  Comparison of 6 machine learning algorithms evaluated on the same test set
                </p>
              </div>

              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b-2 border-gray-300 bg-gray-100">
                      <th className="text-left py-3 px-4 font-semibold">Model</th>
                      <th className="text-center py-3 px-4 font-semibold">Accuracy</th>
                      <th className="text-center py-3 px-4 font-semibold">Precision</th>
                      <th className="text-center py-3 px-4 font-semibold">Recall</th>
                      <th className="text-center py-3 px-4 font-semibold">F1 Score</th>
                      <th className="text-center py-3 px-4 font-semibold">ROC-AUC</th>
                      <th className="text-center py-3 px-4 font-semibold">Training Time</th>
                    </tr>
                  </thead>
                  <tbody>
                    {modelComparisonData.map((model, index) => (
                      <tr
                        key={index}
                        className={`border-b ${
                          model.model === 'XGBoost'
                            ? 'bg-green-50 font-semibold'
                            : 'hover:bg-gray-50'
                        }`}
                      >
                        <td className="py-3 px-4">
                          {model.model}
                          {model.model === 'XGBoost' && (
                            <span className="ml-2 text-xs bg-green-600 text-white px-2 py-1 rounded">
                              SELECTED
                            </span>
                          )}
                        </td>
                        <td className="text-center py-3 px-4">
                          {(model.accuracy * 100).toFixed(1)}%
                        </td>
                        <td className="text-center py-3 px-4">
                          {(model.precision * 100).toFixed(1)}%
                        </td>
                        <td className="text-center py-3 px-4">{(model.recall * 100).toFixed(1)}%</td>
                        <td className="text-center py-3 px-4">{(model.f1 * 100).toFixed(1)}%</td>
                        <td className="text-center py-3 px-4">{(model.rocAuc * 100).toFixed(1)}%</td>
                        <td className="text-center py-3 px-4">{model.trainingTime}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              <div className="bg-blue-50 p-4 rounded-lg border border-blue-200">
                <h5 className="font-semibold text-sm text-gray-800 mb-2">Why XGBoost?</h5>
                <ul className="text-sm text-gray-700 space-y-1 list-disc list-inside">
                  <li>Best overall performance across all metrics</li>
                  <li>Excellent balance between accuracy (83.6%) and speed (2.3s training)</li>
                  <li>Robust to overfitting with built-in regularization</li>
                  <li>Handles imbalanced datasets well (our dataset has slight imbalance)</li>
                  <li>Provides feature importance scores for interpretability</li>
                </ul>
              </div>
            </div>
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  )
}

export default ModelInfo

