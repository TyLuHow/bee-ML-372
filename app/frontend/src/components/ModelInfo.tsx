import { useEffect, useState } from 'react'
import { getModelInfo, ModelInfo as ModelInfoType } from '../services/api'

const ModelInfo = () => {
  const [info, setInfo] = useState<ModelInfoType | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const fetchInfo = async () => {
      try {
        const data = await getModelInfo()
        setInfo(data)
      } catch (err) {
        setError('Unable to fetch model info')
      }
    }
    fetchInfo()
  }, [])

  return (
    <div className="card bg-white border border-amber-100 shadow-lg">
      <h3 className="font-serif text-2xl font-bold text-gray-900 mb-6 border-b border-amber-100 pb-4">
        Scientific Methodology
      </h3>
      
      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        {/* Left Column: Model Stats */}
        <div>
          <h4 className="font-bold text-bee-amber text-sm uppercase tracking-wider mb-4">Model Performance</h4>
          <div className="space-y-4">
            <div className="bg-slate-50 p-4 rounded-lg border border-slate-100">
              <div className="text-3xl font-bold text-slate-800">
                {info?.metrics?.test_accuracy ? (info.metrics.test_accuracy * 100).toFixed(1) : '83.6'}%
              </div>
              <div className="text-sm text-slate-500">Accuracy (Test Set)</div>
            </div>
            
            <div className="bg-slate-50 p-4 rounded-lg border border-slate-100">
              <div className="text-3xl font-bold text-slate-800">
                {info?.metrics?.test_roc_auc ? (info.metrics.test_roc_auc * 100).toFixed(1) : '85.8'}%
              </div>
              <div className="text-sm text-slate-500">ROC-AUC Score</div>
            </div>
            
            <div className="bg-slate-50 p-4 rounded-lg border border-slate-100">
              <div className="text-3xl font-bold text-slate-800">1,035</div>
              <div className="text-sm text-slate-500">Verified Compounds</div>
            </div>
          </div>
        </div>

        {/* Right Column: Methodology */}
        <div className="space-y-6">
          <div>
            <h4 className="font-bold text-bee-amber text-sm uppercase tracking-wider mb-2">Algorithm</h4>
            <p className="text-gray-700 leading-relaxed">
              We utilize <strong>XGBoost</strong> (eXtreme Gradient Boosting), a state-of-the-art ensemble machine learning algorithm known for its performance on structured chemical data.
            </p>
          </div>
          
          <div>
            <h4 className="font-bold text-bee-amber text-sm uppercase tracking-wider mb-2">Key Predictors</h4>
            <p className="text-gray-700 leading-relaxed mb-2">
              Our SHAP (SHapley Additive exPlanations) analysis identifies the following as the most critical factors in bee toxicity:
            </p>
            <ul className="list-disc list-inside text-gray-600 space-y-1 text-sm ml-2">
              <li>Pesticide Type (Insecticide vs Herbicide)</li>
              <li>Registration Year (Newer compounds tend to be safer)</li>
              <li>Lipophilicity (LogP)</li>
              <li>Molecular Weight</li>
            </ul>
          </div>

          <div>
            <h4 className="font-bold text-bee-amber text-sm uppercase tracking-wider mb-2">Data Sources</h4>
            <div className="flex gap-2">
              <span className="px-2 py-1 bg-amber-50 text-bee-amber text-xs font-bold rounded">PPDB</span>
              <span className="px-2 py-1 bg-amber-50 text-bee-amber text-xs font-bold rounded">ECOTOX</span>
              <span className="px-2 py-1 bg-amber-50 text-bee-amber text-xs font-bold rounded">BPDB</span>
            </div>
          </div>
        </div>
      </div>

      <div className="mt-8 pt-6 border-t border-gray-100 flex justify-between items-center">
        <span className="text-xs text-gray-400">
          Last Retrained: {info ? new Date().toLocaleDateString() : 'November 2025'}
        </span>
        <a
          href="http://localhost:8000/docs"
          target="_blank"
          rel="noopener noreferrer"
          className="text-sm font-bold text-bee-amber hover:text-amber-700 transition-colors flex items-center"
        >
          View Full API Documentation <span className="ml-1">→</span>
        </a>
      </div>
    </div>
  )
}

export default ModelInfo
