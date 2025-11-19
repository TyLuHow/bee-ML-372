import { useState } from 'react'
import { predictToxicity, PredictionInput } from '../services/api'
import { Card, CardHeader, CardTitle, CardContent } from './ui/Card'
import { InfoIconTooltip } from './ui/Tooltip'
import { exampleCompounds, exampleOptions } from '../constants/exampleCompounds'
import { validateInput, descriptorTooltips } from '../constants/validationRules'

interface Props {
  onPrediction: (result: any) => void
  onError: (error: string) => void
  onLoading: (loading: boolean) => void
}

interface ValidationErrors {
  [key: string]: string | undefined
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

  const [validationErrors, setValidationErrors] = useState<ValidationErrors>({})
  const [touchedFields, setTouchedFields] = useState<Set<string>>(new Set())
  const [expandedSections, setExpandedSections] = useState<Set<string>>(
    new Set(['compound', 'basic', 'hydrogen', 'structural'])
  )

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()

    // Validate all fields before submission
    const errors: ValidationErrors = {}
    Object.keys(formData).forEach((key) => {
      const value = formData[key as keyof PredictionInput]
      if (typeof value === 'number') {
        const validation = validateInput(key, value)
        if (!validation.isValid) {
          errors[key] = validation.message
        }
      }
    })

    if (Object.keys(errors).length > 0) {
      setValidationErrors(errors)
      onError('Please fix validation errors before submitting.')
      return
    }

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
    const newValue = type === 'number' ? parseFloat(value) : value

    setFormData(prev => ({
      ...prev,
      [name]: newValue
    }))

    // Clear error for this field
    if (validationErrors[name]) {
      setValidationErrors(prev => {
        const newErrors = { ...prev }
        delete newErrors[name]
        return newErrors
      })
    }
  }

  const handleBlur = (e: React.FocusEvent<HTMLInputElement>) => {
    const { name, value } = e.target
    setTouchedFields(prev => new Set(prev).add(name))

    // Validate on blur
    if (value !== '') {
      const numValue = parseFloat(value)
      const validation = validateInput(name, numValue)

      if (!validation.isValid) {
        setValidationErrors(prev => ({
          ...prev,
          [name]: validation.message
        }))
      }
    }
  }

  const handleExampleSelect = (e: React.ChangeEvent<HTMLSelectElement>) => {
    const exampleKey = e.target.value
    if (exampleKey && exampleCompounds[exampleKey]) {
      setFormData(exampleCompounds[exampleKey].data)
      setValidationErrors({})
      setTouchedFields(new Set())
    }
  }

  const toggleSection = (section: string) => {
    setExpandedSections(prev => {
      const newSet = new Set(prev)
      if (newSet.has(section)) {
        newSet.delete(section)
      } else {
        newSet.add(section)
      }
      return newSet
    })
  }

  // Helper component for validated input fields
  const ValidatedInput = ({
    name,
    label,
    tooltip,
    ...props
  }: {
    name: string
    label: string
    tooltip?: string
  } & React.InputHTMLAttributes<HTMLInputElement>) => {
    const error = validationErrors[name]
    const touched = touchedFields.has(name)
    const hasValue = formData[name as keyof PredictionInput] !== undefined
    const isValid = hasValue && !error && touched

    return (
      <div>
        <label className="flex items-center gap-2 text-sm font-medium text-gray-700 mb-1">
          {label}
          {tooltip && <InfoIconTooltip text={tooltip} />}
        </label>
        <input
          name={name}
          value={formData[name as keyof PredictionInput] || ''}
          onChange={handleChange}
          onBlur={handleBlur}
          className={`w-full px-3 py-2 border rounded-lg text-sm transition-colors ${
            error && touched
              ? 'border-red-500 focus:ring-red-500'
              : isValid
              ? 'border-green-500 focus:ring-green-500'
              : 'border-gray-300 focus:ring-indigo-500'
          } focus:outline-none focus:ring-2`}
          {...props}
        />
        {error && touched && (
          <p className="mt-1 text-xs text-red-600">{error}</p>
        )}
        {isValid && (
          <p className="mt-1 text-xs text-green-600 flex items-center gap-1">
            <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
              <path
                fillRule="evenodd"
                d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                clipRule="evenodd"
              />
            </svg>
            Valid
          </p>
        )}
      </div>
    )
  }

  // Helper component for section headers
  const SectionHeader = ({ id, title, icon }: { id: string; title: string; icon: string }) => {
    const isExpanded = expandedSections.has(id)
    return (
      <button
        type="button"
        onClick={() => toggleSection(id)}
        className="w-full flex items-center justify-between p-4 bg-gray-100 hover:bg-gray-200 rounded-lg transition-colors"
      >
        <h3 className="font-semibold text-gray-800 flex items-center gap-2">
          <span>{icon}</span>
          {title}
        </h3>
        <svg
          className={`w-5 h-5 text-gray-600 transition-transform ${
            isExpanded ? 'rotate-180' : ''
          }`}
          fill="none"
          stroke="currentColor"
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>
    )
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle>Enter Compound Properties</CardTitle>
      </CardHeader>
      <CardContent>
        <form onSubmit={handleSubmit} className="space-y-6">
          {/* Load Example Dropdown */}
          <div className="bg-gradient-to-r from-indigo-50 to-purple-50 p-4 rounded-lg border-2 border-indigo-200">
            <label className="block text-sm font-semibold text-gray-800 mb-2">
              Load Example Compound
            </label>
            <select
              onChange={handleExampleSelect}
              className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-indigo-500"
            >
              {exampleOptions.map(option => (
                <option key={option.value} value={option.value}>
                  {option.label}
                </option>
              ))}
            </select>
            <p className="text-xs text-gray-600 mt-2">
              Select a pre-configured example to quickly test the model
            </p>
          </div>

          {/* Section 1: Compound Information */}
          <div className="space-y-3">
            <SectionHeader id="compound" title="Compound Information" icon="📋" />
            {expandedSections.has('compound') && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 p-4 bg-gray-50 rounded-lg">
                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">Source</label>
                  <select
                    name="source"
                    value={formData.source}
                    onChange={handleChange}
                    className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-indigo-500"
                  >
                    <option value="PPDB">PPDB</option>
                    <option value="ECOTOX">ECOTOX</option>
                    <option value="BPDB">BPDB</option>
                  </select>
                </div>

                <ValidatedInput
                  name="year"
                  label="Year"
                  type="number"
                  min={1800}
                  max={2030}
                  tooltip="Year the compound was introduced or studied"
                />

                <div>
                  <label className="block text-sm font-medium text-gray-700 mb-1">
                    Toxicity Type
                  </label>
                  <select
                    name="toxicity_type"
                    value={formData.toxicity_type}
                    onChange={handleChange}
                    className="w-full px-3 py-2 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-indigo-500"
                  >
                    <option value="Contact">Contact</option>
                    <option value="Oral">Oral</option>
                    <option value="Other">Other</option>
                  </select>
                </div>

                <div className="md:col-span-2">
                  <label className="block text-sm font-semibold text-gray-700 mb-2">
                    Chemical Type
                  </label>
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
                    <label className="flex items-center space-x-2 cursor-pointer">
                      <input
                        type="checkbox"
                        name="insecticide"
                        checked={formData.insecticide === 1}
                        onChange={e =>
                          setFormData(prev => ({ ...prev, insecticide: e.target.checked ? 1 : 0 }))
                        }
                        className="w-4 h-4 text-indigo-600"
                      />
                      <span className="text-sm">Insecticide</span>
                    </label>

                    <label className="flex items-center space-x-2 cursor-pointer">
                      <input
                        type="checkbox"
                        name="herbicide"
                        checked={formData.herbicide === 1}
                        onChange={e =>
                          setFormData(prev => ({ ...prev, herbicide: e.target.checked ? 1 : 0 }))
                        }
                        className="w-4 h-4 text-indigo-600"
                      />
                      <span className="text-sm">Herbicide</span>
                    </label>

                    <label className="flex items-center space-x-2 cursor-pointer">
                      <input
                        type="checkbox"
                        name="fungicide"
                        checked={formData.fungicide === 1}
                        onChange={e =>
                          setFormData(prev => ({ ...prev, fungicide: e.target.checked ? 1 : 0 }))
                        }
                        className="w-4 h-4 text-indigo-600"
                      />
                      <span className="text-sm">Fungicide</span>
                    </label>

                    <label className="flex items-center space-x-2 cursor-pointer">
                      <input
                        type="checkbox"
                        name="other_agrochemical"
                        checked={formData.other_agrochemical === 1}
                        onChange={e =>
                          setFormData(prev => ({
                            ...prev,
                            other_agrochemical: e.target.checked ? 1 : 0,
                          }))
                        }
                        className="w-4 h-4 text-indigo-600"
                      />
                      <span className="text-sm">Other</span>
                    </label>
                  </div>
                </div>
              </div>
            )}
          </div>

          {/* Section 2: Basic Molecular Properties */}
          <div className="space-y-3">
            <SectionHeader id="basic" title="Basic Molecular Properties" icon="⚛️" />
            {expandedSections.has('basic') && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 p-4 bg-blue-50 rounded-lg">
                <ValidatedInput
                  name="MolecularWeight"
                  label="Molecular Weight"
                  type="number"
                  step="0.1"
                  tooltip={descriptorTooltips.MolecularWeight}
                />
                <ValidatedInput
                  name="LogP"
                  label="LogP (Lipophilicity)"
                  type="number"
                  step="0.1"
                  tooltip={descriptorTooltips.LogP}
                />
                <ValidatedInput
                  name="TPSA"
                  label="TPSA"
                  type="number"
                  step="0.1"
                  tooltip={descriptorTooltips.TPSA}
                />
                <ValidatedInput
                  name="MolarRefractivity"
                  label="Molar Refractivity"
                  type="number"
                  step="0.1"
                  tooltip={descriptorTooltips.MolarRefractivity}
                />
              </div>
            )}
          </div>

          {/* Section 3: Hydrogen Bonding */}
          <div className="space-y-3">
            <SectionHeader id="hydrogen" title="Hydrogen Bonding" icon="🔗" />
            {expandedSections.has('hydrogen') && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 p-4 bg-green-50 rounded-lg">
                <ValidatedInput
                  name="NumHDonors"
                  label="H-Bond Donors"
                  type="number"
                  tooltip={descriptorTooltips.NumHDonors}
                />
                <ValidatedInput
                  name="NumHAcceptors"
                  label="H-Bond Acceptors"
                  type="number"
                  tooltip={descriptorTooltips.NumHAcceptors}
                />
              </div>
            )}
          </div>

          {/* Section 4: Structural Features */}
          <div className="space-y-3">
            <SectionHeader id="structural" title="Structural Features" icon="🧬" />
            {expandedSections.has('structural') && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 p-4 bg-purple-50 rounded-lg">
                <ValidatedInput
                  name="NumRotatableBonds"
                  label="Rotatable Bonds"
                  type="number"
                  tooltip={descriptorTooltips.NumRotatableBonds}
                />
                <ValidatedInput
                  name="AromaticRings"
                  label="Aromatic Rings"
                  type="number"
                  tooltip={descriptorTooltips.AromaticRings}
                />
                <ValidatedInput
                  name="NumAromaticRings"
                  label="Num Aromatic Rings"
                  type="number"
                  tooltip={descriptorTooltips.NumAromaticRings}
                />
                <ValidatedInput
                  name="NumRings"
                  label="Total Rings"
                  type="number"
                  tooltip={descriptorTooltips.NumRings}
                />
                <ValidatedInput
                  name="RingCount"
                  label="Ring Count"
                  type="number"
                  tooltip={descriptorTooltips.RingCount}
                />
                <ValidatedInput
                  name="NumSaturatedRings"
                  label="Saturated Rings"
                  type="number"
                  tooltip={descriptorTooltips.NumSaturatedRings}
                />
                <ValidatedInput
                  name="NumAliphaticRings"
                  label="Aliphatic Rings"
                  type="number"
                  tooltip={descriptorTooltips.NumAliphaticRings}
                />
                <ValidatedInput
                  name="NumAromaticCarbocycles"
                  label="Aromatic Carbocycles"
                  type="number"
                  tooltip={descriptorTooltips.NumAromaticCarbocycles}
                />
                <ValidatedInput
                  name="NumSaturatedCarbocycles"
                  label="Saturated Carbocycles"
                  type="number"
                  tooltip={descriptorTooltips.NumSaturatedCarbocycles}
                />
                <ValidatedInput
                  name="NumHeteroatoms"
                  label="Heteroatoms"
                  type="number"
                  tooltip={descriptorTooltips.NumHeteroatoms}
                />
                <ValidatedInput
                  name="NumAromaticAtoms"
                  label="Aromatic Atoms"
                  type="number"
                  tooltip={descriptorTooltips.NumAromaticAtoms}
                />
                <ValidatedInput
                  name="HeavyAtomCount"
                  label="Heavy Atoms"
                  type="number"
                  tooltip={descriptorTooltips.HeavyAtomCount}
                />
                <ValidatedInput
                  name="FractionCsp3"
                  label="Fraction Csp3"
                  type="number"
                  step="0.01"
                  tooltip={descriptorTooltips.FractionCsp3}
                />
                <ValidatedInput
                  name="FractionCSP3"
                  label="Fraction CSP3"
                  type="number"
                  step="0.01"
                  tooltip={descriptorTooltips.FractionCSP3}
                />
                <ValidatedInput
                  name="BertzCT"
                  label="Bertz Complexity"
                  type="number"
                  step="0.1"
                  tooltip={descriptorTooltips.BertzCT}
                />
              </div>
            )}
          </div>

          {/* Submit Button */}
          <button
            type="submit"
            className="w-full bg-gradient-to-r from-indigo-600 to-purple-600 hover:from-indigo-700 hover:to-purple-700 text-white font-semibold py-3 px-6 rounded-lg shadow-lg transition-all duration-200 transform hover:scale-105"
          >
            Predict Toxicity
          </button>

          <p className="text-xs text-gray-500 text-center">
            Tip: Insecticides are typically more toxic to bees than herbicides or fungicides
          </p>
        </form>
      </CardContent>
    </Card>
  )
}

export default PredictionForm

