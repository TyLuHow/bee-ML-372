import React, { useState } from 'react'
import PredictionForm from '../components/PredictionForm'
import ResultDisplay from '../components/ResultDisplay'
import { Card, CardContent } from '../components/ui/Card'
import { useAppContext, PredictionResult } from '../store/AppContext'

export const PredictPage: React.FC = () => {
  const { addPrediction } = useAppContext()
  const [result, setResult] = useState<PredictionResult | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const handlePrediction = (predictionResult: PredictionResult) => {
    setResult(predictionResult)
    setError(null)
    addPrediction(predictionResult)
  }

  const handleError = (errorMessage: string) => {
    setError(errorMessage)
    setResult(null)
  }

  const handleLoading = (isLoading: boolean) => {
    setLoading(isLoading)
  }

  return (
    <div className="space-y-6">
      {/* Instructions */}
      <Card>
        <CardContent className="py-4">
          <div className="flex items-start gap-3">
            <span className="text-2xl">💡</span>
            <div>
              <p className="text-text-primary font-medium mb-1">How to use</p>
              <p className="text-sm text-text-secondary">
                Enter the molecular properties of a pesticide compound below. The ML model will predict
                whether it's toxic or non-toxic to honey bees.
              </p>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Main Content Grid */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Prediction Form - Takes 2 columns */}
        <div className="lg:col-span-2">
          <PredictionForm
            onPrediction={handlePrediction}
            onError={handleError}
            onLoading={handleLoading}
          />
        </div>

        {/* Results Display - Takes 1 column */}
        <div>
          <ResultDisplay
            result={result}
            loading={loading}
            error={error}
          />
        </div>
      </div>
    </div>
  )
}
