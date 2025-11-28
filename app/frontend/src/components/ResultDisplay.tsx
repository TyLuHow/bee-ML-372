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

const ResultDisplay = ({ result, loading, error }: Props) => {
  if (loading) {
    return (
      <div className="card bg-white/80 backdrop-blur text-center py-16 border-2 border-amber-100">
        <div className="relative mx-auto h-20 w-20">
          <div className="absolute inset-0 border-t-4 border-bee-amber rounded-full animate-spin"></div>
          <div className="absolute inset-0 flex items-center justify-center text-2xl">🐝</div>
        </div>
        <h3 className="mt-6 text-xl font-serif font-bold text-gray-800">Processing Compound</h3>
        <p className="text-gray-500 mt-2">Analyzing molecular structure against 1,000+ known toxins...</p>
      </div>
    )
  }

  if (error) {
    return (
      <div className="card bg-red-50 border-l-4 border-red-500 p-6 shadow-lg">
        <div className="flex items-start">
          <span className="text-3xl mr-4">⚠️</span>
          <div>
            <h3 className="text-lg font-bold text-red-900 mb-1">Analysis Failed</h3>
            <p className="text-red-700 text-sm leading-relaxed">{error}</p>
            <p className="text-xs text-red-500 mt-3 font-mono bg-red-100 p-2 rounded inline-block">
              Check API connection at port 8001
            </p>
          </div>
        </div>
      </div>
    )
  }

  if (!result) {
    return (
      <div className="card bg-gradient-to-br from-white to-amber-50/50 border-2 border-dashed border-amber-200 text-center py-16 px-6">
        <div className="text-6xl mb-6 opacity-50 grayscale hover:grayscale-0 transition-all duration-500 transform hover:scale-110 cursor-default">
          🐝
        </div>
        <h3 className="text-xl font-serif font-bold text-gray-800 mb-2">Awaiting Input</h3>
        <p className="text-gray-500 text-sm max-w-xs mx-auto">
          The bees are waiting. Enter chemical properties to assess potential toxicity risks.
        </p>
      </div>
    )
  }

  const isToxic = result.prediction === 1
  const confidencePct = Math.round(result.confidence * 100)
  
  // Calculate probabilities with fallbacks
  const probSafe = result.probability_non_toxic ?? (1 - result.confidence)
  const probToxic = result.probability_toxic ?? result.confidence
  const safePercent = Math.max(0, Math.min(100, probSafe * 100))
  const toxicPercent = Math.max(0, Math.min(100, probToxic * 100))
  
  // Theme configuration based on result
  const theme = isToxic 
    ? {
        bg: 'bg-rose-50',
        border: 'border-rose-200',
        text: 'text-rose-900',
        accent: 'bg-rose-500',
        icon: '⚠️',
        title: 'High Toxicity Risk',
        desc: 'This compound shares structural alerts with known bee toxins.'
      }
    : {
        bg: 'bg-emerald-50',
        border: 'border-emerald-200',
        text: 'text-emerald-900',
        accent: 'bg-emerald-500',
        icon: '🌿',
        title: 'Likely Safe',
        desc: 'Molecular profile suggests low risk to honey bees.'
      }

  return (
    <div className={`card ${theme.bg} border-2 ${theme.border} shadow-xl overflow-hidden relative`}>
      {/* Result Header */}
      <div className="text-center pb-8 border-b border-black/5 relative z-10">
        <div className="text-6xl mb-4 filter drop-shadow-sm">{theme.icon}</div>
        <h2 className={`text-3xl font-serif font-bold ${theme.text} mb-2`}>
          {theme.title}
        </h2>
        <p className={`${theme.text} opacity-80 text-sm max-w-xs mx-auto leading-relaxed`}>
          {theme.desc}
        </p>
      </div>

      {/* Confidence Meter */}
      <div className="py-8 px-4 bg-white/50">
        <div className="flex justify-between items-end mb-2">
          <span className="text-xs font-bold uppercase tracking-widest text-gray-500">
            AI Confidence
          </span>
          <span className="text-2xl font-bold text-gray-800">
            {confidencePct}%
          </span>
        </div>
        
        {/* Probability Bar */}
        <div className="h-5 bg-gray-100 rounded-full overflow-hidden flex shadow-inner">
          <div 
            className="h-full transition-all duration-1000 ease-out relative group"
            style={{ 
              width: `${safePercent}%`,
              background: 'linear-gradient(90deg, #10b981 0%, #34d399 100%)'
            }}
          >
            {safePercent > 15 && (
              <span className="absolute inset-0 flex items-center justify-center text-white text-xs font-bold drop-shadow">
                {safePercent.toFixed(0)}%
              </span>
            )}
          </div>
          <div 
            className="h-full transition-all duration-1000 ease-out relative group"
            style={{ 
              width: `${toxicPercent}%`,
              background: 'linear-gradient(90deg, #f87171 0%, #ef4444 100%)'
            }}
          >
            {toxicPercent > 15 && (
              <span className="absolute inset-0 flex items-center justify-center text-white text-xs font-bold drop-shadow">
                {toxicPercent.toFixed(0)}%
              </span>
            )}
          </div>
        </div>
        <div className="flex justify-between mt-2 text-xs text-gray-400 font-medium">
          <span>Safe</span>
          <span>Toxic</span>
        </div>
      </div>

      {/* Action Item / Context */}
      <div className="bg-white/80 p-5 text-sm">
        <h4 className="font-bold text-gray-700 mb-2 flex items-center">
          <span className="mr-2 text-lg">💡</span>
          Recommendation
        </h4>
        <p className="text-gray-600 leading-relaxed">
          {isToxic 
            ? "Avoid use during bloom. Consider searching for alternative compounds with lower toxicity profiles in the BPDB (Bio-Pesticide) database."
            : "Proceed with standard field testing protocols. While AI prediction is favorable, always verify with in vivo assays for final regulatory approval."
          }
        </p>
      </div>

      {/* Footer */}
      <div className="bg-gray-50 py-3 px-5 text-[10px] text-gray-400 flex justify-between items-center">
        <span>Model: XGBoost v1.0</span>
        <span>{new Date(result.timestamp).toLocaleTimeString()}</span>
      </div>
    </div>
  )
}

export default ResultDisplay
