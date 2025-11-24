import { useState } from 'react'
import PredictionForm from './components/PredictionForm'
import ResultDisplay from './components/ResultDisplay'
import ModelInfo from './components/ModelInfo'
import './App.css'

interface PredictionResult {
  prediction: number
  label_text: string
  probability_toxic: number
  probability_non_toxic: number
  confidence: number
  timestamp: string
}

function App() {
  const [result, setResult] = useState<PredictionResult | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [activeTab, setActiveTab] = useState<'predict' | 'model'>('predict')

  const handlePrediction = (predictionResult: PredictionResult) => {
    setResult(predictionResult)
    setError(null)
    // Scroll to result on mobile
    const resultElement = document.getElementById('result-section')
    if (resultElement) {
      resultElement.scrollIntoView({ behavior: 'smooth' })
    }
  }

  const handleError = (errorMessage: string) => {
    setError(errorMessage)
    setResult(null)
  }

  return (
    <div className="min-h-screen bg-nature-cream bg-honeycomb">
      {/* Hero Section */}
      <div className="bg-gradient-to-b from-white via-amber-50/50 to-transparent pt-16 pb-12 border-b border-amber-100/50">
        <div className="max-w-6xl mx-auto px-4 sm:px-6 text-center">
          <div className="inline-flex items-center justify-center p-2 bg-amber-100 rounded-full mb-6 animate-pulse">
            <span className="text-2xl mr-2">🐝</span>
            <span className="text-bee-amber font-medium text-sm pr-2">Project ApisTox</span>
          </div>
          <h1 className="text-5xl md:text-7xl font-bold text-bee-black mb-6 tracking-tight leading-tight">
            Protecting Pollinators <br/>
            <span className="text-transparent bg-clip-text bg-gradient-to-r from-bee-gold to-bee-amber">
              Before It's Too Late
            </span>
          </h1>
          <p className="text-xl text-gray-600 max-w-2xl mx-auto leading-relaxed font-light">
            Pesticides are a leading cause of colony collapse. We use advanced machine learning to predict chemical toxicity to honey bees—empowering safer choices for agriculture and our planet.
          </p>
        </div>
      </div>

      {/* Navigation Tabs */}
      <div className="max-w-6xl mx-auto px-4 mt-8 mb-8">
        <div className="flex justify-center space-x-4 border-b border-amber-200/50 pb-1">
          <button 
            onClick={() => setActiveTab('predict')}
            className={`pb-3 px-4 text-lg font-medium transition-colors relative ${
              activeTab === 'predict' 
                ? 'text-bee-amber' 
                : 'text-gray-500 hover:text-gray-700'
            }`}
          >
            Prediction Engine
            {activeTab === 'predict' && (
              <span className="absolute bottom-0 left-0 w-full h-0.5 bg-bee-amber rounded-t-full"></span>
            )}
          </button>
          <button 
            onClick={() => setActiveTab('model')}
            className={`pb-3 px-4 text-lg font-medium transition-colors relative ${
              activeTab === 'model' 
                ? 'text-bee-amber' 
                : 'text-gray-500 hover:text-gray-700'
            }`}
          >
            Scientific Basis
            {activeTab === 'model' && (
              <span className="absolute bottom-0 left-0 w-full h-0.5 bg-bee-amber rounded-t-full"></span>
            )}
          </button>
        </div>
      </div>

      {/* Main Application Area */}
      <main className="max-w-7xl mx-auto px-4 pb-24">
        {activeTab === 'predict' ? (
          <div className="grid grid-cols-1 lg:grid-cols-12 gap-8 items-start">
            {/* Input Column */}
            <div className="lg:col-span-7 xl:col-span-8">
              <PredictionForm
                onPrediction={handlePrediction}
                onError={handleError}
                onLoading={setLoading}
              />
            </div>

            {/* Result Column */}
            <div className="lg:col-span-5 xl:col-span-4 sticky top-8" id="result-section">
              <ResultDisplay
                result={result}
                loading={loading}
                error={error}
              />
              
              {/* Quick Context */}
              {!loading && !result && (
                <div className="mt-6 bg-white/50 p-6 rounded-2xl border border-amber-100 backdrop-blur-sm">
                  <h3 className="font-serif text-lg font-bold text-gray-800 mb-3">Why prediction matters</h3>
                  <p className="text-gray-600 text-sm leading-relaxed">
                    Traditional toxicity testing is slow and expensive. In silico (computational) models allow us to screen thousands of compounds instantly, flagging potential risks before they ever reach the field.
                  </p>
                </div>
              )}
            </div>
          </div>
        ) : (
          <div className="max-w-4xl mx-auto">
            <ModelInfo />
          </div>
        )}
      </main>

      {/* Footer */}
      <footer className="bg-white border-t border-amber-100 mt-12 py-12">
        <div className="max-w-6xl mx-auto px-4 text-center">
          <p className="font-serif text-2xl text-bee-black mb-4">ApisTox</p>
          <div className="flex justify-center space-x-8 mb-8 text-sm text-gray-500">
            <span>XGBoost Classifier</span>
            <span>•</span>
            <span>1,035 Compounds Analyzed</span>
            <span>•</span>
            <span>SHAP Interpretability</span>
          </div>
          <p className="text-gray-400 text-sm">
            Built for IME 372 • Cal Poly 2025
          </p>
        </div>
      </footer>
    </div>
  )
}

export default App
