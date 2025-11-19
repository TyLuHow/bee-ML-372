import React from 'react'
import { Card, CardHeader, CardTitle, CardContent } from '../components/ui/Card'
import { Badge } from '../components/ui/Badge'
import { useAppContext } from '../store/AppContext'

export const Dashboard: React.FC = () => {
  const { predictions } = useAppContext()

  const recentPredictions = predictions.slice(0, 5)
  const toxicCount = predictions.filter(p => p.prediction === 1).length
  const safeCount = predictions.filter(p => p.prediction === 0).length

  return (
    <div className="space-y-6">
      {/* Stats Grid */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <Card>
          <CardContent className="pt-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-text-secondary">Total Predictions</p>
                <p className="text-3xl font-bold text-text-primary mt-1">{predictions.length}</p>
              </div>
              <div className="text-4xl">🎯</div>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardContent className="pt-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-text-secondary">Toxic Compounds</p>
                <p className="text-3xl font-bold text-accent-toxic mt-1">{toxicCount}</p>
              </div>
              <div className="text-4xl">⚠️</div>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardContent className="pt-6">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-text-secondary">Safe Compounds</p>
                <p className="text-3xl font-bold text-accent-safe mt-1">{safeCount}</p>
              </div>
              <div className="text-4xl">✅</div>
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Recent Predictions */}
      <Card>
        <CardHeader>
          <CardTitle>Recent Predictions</CardTitle>
        </CardHeader>
        <CardContent>
          {recentPredictions.length === 0 ? (
            <div className="text-center py-8 text-text-secondary">
              <p className="text-4xl mb-2">🐝</p>
              <p>No predictions yet. Navigate to Predict Toxicity to get started!</p>
            </div>
          ) : (
            <div className="space-y-3">
              {recentPredictions.map((pred, index) => (
                <div
                  key={index}
                  className="flex items-center justify-between p-4 bg-bg-secondary rounded-lg"
                >
                  <div className="flex-1">
                    <p className="font-medium text-text-primary">
                      {pred.compound_name || `Prediction ${index + 1}`}
                    </p>
                    <p className="text-sm text-text-tertiary">
                      {new Date(pred.timestamp).toLocaleString()}
                    </p>
                  </div>
                  <div className="flex items-center gap-3">
                    <Badge variant={pred.prediction === 1 ? 'toxic' : 'safe'}>
                      {pred.label_text}
                    </Badge>
                    <span className="text-sm text-text-secondary">
                      {(pred.confidence * 100).toFixed(1)}% confidence
                    </span>
                  </div>
                </div>
              ))}
            </div>
          )}
        </CardContent>
      </Card>

      {/* Quick Info */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <Card>
          <CardHeader>
            <CardTitle>About ApisTox</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-text-secondary mb-4">
              This application uses machine learning to predict the toxicity of pesticides to honey bees,
              helping protect our vital pollinators.
            </p>
            <div className="space-y-2 text-sm">
              <div className="flex justify-between">
                <span className="text-text-secondary">Dataset:</span>
                <span className="font-medium">1,035 compounds</span>
              </div>
              <div className="flex justify-between">
                <span className="text-text-secondary">Model:</span>
                <span className="font-medium">XGBoost Classifier</span>
              </div>
              <div className="flex justify-between">
                <span className="text-text-secondary">Accuracy:</span>
                <span className="font-medium text-accent-safe">83.6%</span>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Quick Links</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              <a
                href="#/predict"
                className="flex items-center gap-3 p-3 bg-bg-secondary rounded-lg hover:bg-bg-tertiary transition-colors"
              >
                <span className="text-2xl">🎯</span>
                <div>
                  <p className="font-medium text-text-primary">Make a Prediction</p>
                  <p className="text-sm text-text-tertiary">Test a compound for toxicity</p>
                </div>
              </a>
              <a
                href="#/explorer"
                className="flex items-center gap-3 p-3 bg-bg-secondary rounded-lg hover:bg-bg-tertiary transition-colors"
              >
                <span className="text-2xl">🔍</span>
                <div>
                  <p className="font-medium text-text-primary">Explore Dataset</p>
                  <p className="text-sm text-text-tertiary">Browse the ApisTox data</p>
                </div>
              </a>
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
