import { useState } from 'react'
import { predictToxicity, PredictionInput } from '../services/api'

interface Props {
  onPrediction: (result: any) => void
  onError: (error: string) => void
  onLoading: (loading: boolean) => void
}

const PredictionForm = ({ onPrediction, onError, onLoading }: Props) => {
  const [formData, setFormData] = useState<Partial<PredictionInput>>({
    source: 'PPDB',
    year: 2020,
    toxicity_type: 'Contact',
    insecticide: 1,
    herbicide: 0,
    fungicide: 0,
    other_agrochemical: 0,
    MolecularWeight: 350.5,
    LogP: 3.2,
    NumHDonors: 2,
    NumHAcceptors: 4,
    NumRotatableBonds: 5,
    AromaticRings: 2,
    NumAromaticRings: 2,
    TPSA: 65.3,
    NumHeteroatoms: 5,
    NumAromaticAtoms: 12,
    NumSaturatedRings: 0,
    NumAliphaticRings: 0,
    RingCount: 2,
    NumRings: 2,
    FractionCsp3: 0.25,
    FractionCSP3: 0.25,
    NumAromaticCarbocycles: 1,
    NumSaturatedCarbocycles: 0,
    MolarRefractivity: 95.5,
    BertzCT: 850.0,
    HeavyAtomCount: 25
  })

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    onLoading(true)
    try {
      const result = await predictToxicity(formData as PredictionInput)
      onPrediction(result)
    } catch (err: any) {
      onError(err.response?.data?.detail || 'Prediction failed. Please check your inputs.')
    } finally {
      onLoading(false)
    }
  }

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
    const { name, value, type } = e.target
    setFormData(prev => ({
      ...prev,
      [name]: type === 'number' ? parseFloat(value) : value
    }))
  }

  const toggleCheckbox = (name: string) => {
    setFormData(prev => ({
      ...prev,
      // @ts-ignore
      [name]: prev[name] === 1 ? 0 : 1
    }))
  }

  return (
    <div className="card bg-white/90 backdrop-blur border-amber-100/50 shadow-xl shadow-amber-100/20">
      <div className="mb-8">
        <h2 className="text-3xl font-serif font-bold text-gray-900 mb-2">
          Analyze a Compound
        </h2>
        <p className="text-gray-500">
          Input the chemical properties below to generate a real-time toxicity risk assessment.
        </p>
      </div>
      
      <form onSubmit={handleSubmit} className="space-y-8">
        
        {/* Section 1: Identity & Classification */}
        <section>
          <h3 className="text-sm uppercase tracking-wider text-gray-500 font-bold mb-4 border-b border-gray-100 pb-2">
            1. Chemical Identity
          </h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
            <div>
              <label className="label text-gray-600">Data Source</label>
              <select name="source" value={formData.source} onChange={handleChange} className="input-field bg-gray-50/50">
                <option value="PPDB">PPDB (Pesticide Properties)</option>
                <option value="ECOTOX">ECOTOX Database</option>
                <option value="BPDB">Bio-Pesticide DB</option>
              </select>
            </div>
            
            <div>
              <label className="label text-gray-600">Registration Year</label>
              <input
                type="number"
                name="year"
                value={formData.year}
                onChange={handleChange}
                className="input-field bg-gray-50/50"
                min="1800"
                max="2030"
              />
            </div>
            
            <div>
              <label className="label text-gray-600">Exposure Route</label>
              <select name="toxicity_type" value={formData.toxicity_type} onChange={handleChange} className="input-field bg-gray-50/50">
                <option value="Contact">Contact (Direct Spray)</option>
                <option value="Oral">Oral (Ingestion)</option>
                <option value="Other">Other</option>
              </select>
            </div>
          </div>
        </section>

        {/* Section 2: Usage Type */}
        <section>
          <h3 className="text-sm uppercase tracking-wider text-gray-500 font-bold mb-4 border-b border-gray-100 pb-2">
            2. Pesticide Class
          </h3>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            {[
              { key: 'insecticide', label: 'Insecticide', icon: '🦗' },
              { key: 'herbicide', label: 'Herbicide', icon: '🌱' },
              { key: 'fungicide', label: 'Fungicide', icon: '🍄' },
              { key: 'other_agrochemical', label: 'Other', icon: '🧪' }
            ].map((type) => (
              <div 
                key={type.key}
                // @ts-ignore
                onClick={() => toggleCheckbox(type.key)}
                // @ts-ignore
                className={`cursor-pointer rounded-xl border-2 p-4 flex flex-col items-center justify-center transition-all duration-200 ${formData[type.key] === 1 ? 'border-bee-gold bg-amber-50 text-bee-amber' : 'border-gray-100 hover:border-gray-200 text-gray-500'}`}
              >
                <span className="text-2xl mb-2">{type.icon}</span>
                <span className="font-medium text-sm">{type.label}</span>
              </div>
            ))}
          </div>
        </section>

        {/* Section 3: Molecular Properties */}
        <section>
          <div className="flex justify-between items-end mb-4 border-b border-gray-100 pb-2">
            <h3 className="text-sm uppercase tracking-wider text-gray-500 font-bold">
              3. Molecular Fingerprint
            </h3>
            <span className="text-xs text-bee-amber bg-amber-50 px-2 py-1 rounded-full">
              Crucial for prediction
            </span>
          </div>
          
          <div className="bg-slate-50/50 p-6 rounded-2xl border border-slate-100">
            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-x-6 gap-y-6">
              
              {/* Key properties first */}
              <div className="col-span-2 md:col-span-1">
                <label className="label text-bee-black font-semibold">Molecular Weight</label>
                <div className="relative">
                  <input type="number" step="0.1" name="MolecularWeight" value={formData.MolecularWeight} onChange={handleChange} className="input-field pl-3" />
                  <span className="absolute right-3 top-2 text-gray-400 text-xs">g/mol</span>
                </div>
                <p className="text-[10px] text-gray-400 mt-1">Heavier molecules may have lower bioavailability.</p>
              </div>
              
              <div className="col-span-2 md:col-span-1">
                <label className="label text-bee-black font-semibold">LogP (Lipophilicity)</label>
                <input type="number" step="0.1" name="LogP" value={formData.LogP} onChange={handleChange} className="input-field" />
                <p className="text-[10px] text-gray-400 mt-1">High LogP often means higher toxicity accumulation.</p>
              </div>

              <div className="col-span-2 md:col-span-1">
                <label className="label text-bee-black font-semibold">TPSA</label>
                <input type="number" step="0.1" name="TPSA" value={formData.TPSA} onChange={handleChange} className="input-field" />
                <p className="text-[10px] text-gray-400 mt-1">Polar surface area affects transport.</p>
              </div>

              {/* Secondary properties */}
              <div className="border-t border-slate-200 col-span-full my-2 pt-2"></div>

              <div>
                <label className="label text-xs">H-Bond Donors</label>
                <input type="number" name="NumHDonors" value={formData.NumHDonors} onChange={handleChange} className="input-field text-sm py-1" />
              </div>
              
              <div>
                <label className="label text-xs">H-Bond Acceptors</label>
                <input type="number" name="NumHAcceptors" value={formData.NumHAcceptors} onChange={handleChange} className="input-field text-sm py-1" />
              </div>
              
              <div>
                <label className="label text-xs">Rotatable Bonds</label>
                <input type="number" name="NumRotatableBonds" value={formData.NumRotatableBonds} onChange={handleChange} className="input-field text-sm py-1" />
              </div>
              
              <div>
                <label className="label text-xs">Aromatic Rings</label>
                <input type="number" name="AromaticRings" value={formData.AromaticRings} onChange={handleChange} className="input-field text-sm py-1" />
              </div>
              
              <div>
                <label className="label text-xs">Heteroatoms</label>
                <input type="number" name="NumHeteroatoms" value={formData.NumHeteroatoms} onChange={handleChange} className="input-field text-sm py-1" />
              </div>
              
              <div>
                <label className="label text-xs">Heavy Atoms</label>
                <input type="number" name="HeavyAtomCount" value={formData.HeavyAtomCount} onChange={handleChange} className="input-field text-sm py-1" />
              </div>
            </div>
          </div>
        </section>

        <button
          type="submit"
          className="w-full bg-gradient-to-r from-bee-gold to-bee-amber hover:from-bee-amber hover:to-amber-700 text-white font-serif font-bold text-xl py-4 rounded-xl shadow-lg shadow-amber-500/30 transform hover:-translate-y-1 transition-all duration-200 flex items-center justify-center gap-3"
        >
          <span>Run Predictive Analysis</span>
          <span className="text-2xl">🔮</span>
        </button>
      </form>
    </div>
  )
}

export default PredictionForm
