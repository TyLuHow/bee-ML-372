import React, { useState } from 'react';
import {
  Leaf,
  ArrowRight,
  Activity,
  AlertTriangle,
  CheckCircle2,
  Info,
  Beaker,
  BookOpen,
  Microscope,
  Flower2,
  Sprout,
  Database,
  Search,
  ChevronRight,
  ChevronDown,
  Play,
  BarChart3,
  FileText,
  Atom
} from 'lucide-react';
import { ChemicalData, PredictionResult, TabView, Scenario, CompoundData } from './types';
import { analyzeChemicalToxicity } from './services/geminiService';

// --- Shared Components ---

const Button = ({ children, onClick, variant = 'primary', className = '', disabled = false }: any) => {
  const base = "px-6 py-3 rounded-xl font-serif font-bold transition-all duration-quick disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center gap-2";
  const styles = {
    primary: "bg-journal-accent text-white shadow-card-subtle hover:shadow-lifted hover:bg-[#B45309]",
    secondary: "bg-white border border-journal-border text-journal-text hover:border-journal-accent hover:text-journal-accent shadow-paper",
    ghost: "bg-transparent text-journal-dim hover:text-journal-accent"
  };
  return (
    <button onClick={onClick} disabled={disabled} className={`${base} ${styles[variant as keyof typeof styles]} ${className}`}>
      {children}
    </button>
  );
};

const Card = ({ children, className = '' }: any) => (
  <div className={`bg-white border border-[#D4CFC5] rounded-xl shadow-card-subtle p-8 transition-all duration-smooth ${className}`}>
    {children}
  </div>
);

const SectionHeader = ({ title, subtitle }: { title: string, subtitle?: string }) => (
  <div className="mb-10 border-l-4 border-journal-accent pl-6">
    <h2 className="text-2xl font-display font-bold text-journal-text">{title}</h2>
    {subtitle && <p className="text-journal-dim font-serif italic mt-2">{subtitle}</p>}
  </div>
);

// --- Page Views ---

const PredictionEngine = () => {
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<PredictionResult | null>(null);
  const [expandedSection, setExpandedSection] = useState<'basic' | 'descriptors' | null>('descriptors');
  const [formData, setFormData] = useState<ChemicalData>({
    name: '',
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
    NumAromaticRings: 2,
    TPSA: 65.3,
    NumHeteroatoms: 5,
    NumRings: 2,
    NumSaturatedRings: 0,
    NumAliphaticRings: 0,
    FractionCSP3: 0.25,
    MolarRefractivity: 95.5,
    BertzCT: 850.0,
    HeavyAtomCount: 25
  });

  const handleAnalyze = async () => {
    setLoading(true);
    const res = await analyzeChemicalToxicity(formData);
    setResult(res);
    setLoading(false);
  };

  const updateField = (field: keyof ChemicalData, value: any) => {
    setFormData(prev => ({ ...prev, [field]: value }));
  };

  const setCategory = (cat: string) => {
    updateField('category', cat);
    updateField('insecticide', cat === 'Insecticide' ? 1 : 0);
    updateField('herbicide', cat === 'Herbicide' ? 1 : 0);
    updateField('fungicide', cat === 'Fungicide' ? 1 : 0);
    updateField('other_agrochemical', cat === 'Other' ? 1 : 0);
  };

  return (
    <div className="animate-fade-in">
      <SectionHeader
        title="Analyze a Compound"
        subtitle="Input chemical properties below to generate a real-time toxicity risk assessment."
      />

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-8 items-start">
        {/* Input Form */}
        <div className="lg:col-span-7">
          <Card className="space-y-6">
            {/* Basic Information */}
            <div>
              <button
                onClick={() => setExpandedSection(expandedSection === 'basic' ? null : 'basic')}
                className="w-full flex justify-between items-center border-b border-journal-border pb-4 mb-4 hover:text-journal-accent transition-colors"
              >
                <h3 className="font-bold text-xs uppercase tracking-widest text-journal-dim">1. Chemical Identity & Classification</h3>
                {expandedSection === 'basic' ? <ChevronDown size={20} /> : <ChevronRight size={20} />}
              </button>

              {expandedSection === 'basic' && (
                <div className="space-y-4 animate-fade-in">
                  <div>
                    <label className="block text-sm font-bold text-journal-text mb-2">Compound Name</label>
                    <div className="relative">
                      <input
                        type="text"
                        placeholder="e.g. Imidacloprid"
                        value={formData.name || ''}
                        onChange={(e) => updateField('name', e.target.value)}
                        className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg focus:border-journal-accent focus:bg-white focus:shadow-sm outline-none transition-all duration-quick font-medium"
                      />
                      <Search size={18} className="absolute right-4 top-3.5 text-journal-dim/40" />
                    </div>
                  </div>

                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <label className="block text-sm font-bold text-journal-text mb-2">Data Source</label>
                      <select
                        value={formData.source || 'PPDB'}
                        onChange={(e) => updateField('source', e.target.value)}
                        className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg outline-none text-sm font-medium focus:border-journal-accent focus:bg-white focus:shadow-sm transition-all duration-quick cursor-pointer"
                      >
                        <option value="PPDB">PPDB (Pesticide Properties)</option>
                        <option value="ECOTOX">ECOTOX (EPA)</option>
                        <option value="BPDB">BPDB (Bio-Pesticides)</option>
                        <option value="Unknown">Unknown</option>
                      </select>
                    </div>
                    <div>
                      <label className="block text-sm font-bold text-journal-text mb-2">Registration Year</label>
                      <input
                        type="number"
                        value={formData.year || 2020}
                        onChange={(e) => updateField('year', parseInt(e.target.value))}
                        className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg focus:border-journal-accent focus:bg-white focus:shadow-sm outline-none transition-all duration-quick font-mono"
                        min="1800"
                        max="2030"
                      />
                    </div>
                  </div>

                  <div>
                    <label className="block text-sm font-bold text-journal-text mb-2">Exposure Route / Toxicity Type</label>
                    <select
                      value={formData.toxicity_type || 'Contact'}
                      onChange={(e) => updateField('toxicity_type', e.target.value)}
                      className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg outline-none text-sm font-medium focus:border-journal-accent focus:bg-white focus:shadow-sm transition-all duration-quick cursor-pointer"
                    >
                      <option value="Contact">Contact (Direct Spray)</option>
                      <option value="Oral">Oral (Ingestion)</option>
                      <option value="Systemic">Systemic</option>
                      <option value="Other">Other</option>
                    </select>
                  </div>

                  <div>
                    <label className="block text-sm font-bold text-journal-text mb-3">Pesticide Class</label>
                    <div className="grid grid-cols-2 gap-3">
                      {['Insecticide', 'Herbicide', 'Fungicide', 'Other'].map(cat => (
                        <button
                          key={cat}
                          onClick={() => setCategory(cat)}
                          className={`p-5 border-2 rounded-xl flex flex-col items-center gap-2 transition-all duration-quick ${
                            formData.category === cat ||
                            (cat === 'Insecticide' && formData.insecticide === 1) ||
                            (cat === 'Herbicide' && formData.herbicide === 1) ||
                            (cat === 'Fungicide' && formData.fungicide === 1) ||
                            (cat === 'Other' && formData.other_agrochemical === 1)
                            ? 'border-journal-accent bg-[#FEF3C7] text-journal-accent shadow-card-subtle scale-[1.02]'
                            : 'border-journal-border hover:border-journal-accent/60 hover:bg-[#FFFAF0]'
                          }`}
                        >
                          {cat === 'Insecticide' && <BugIcon />}
                          {cat === 'Herbicide' && <Sprout size={22} />}
                          {cat === 'Fungicide' && <Leaf size={22} />}
                          {cat === 'Other' && <Beaker size={22} />}
                          <span className="text-sm font-bold">{cat}</span>
                        </button>
                      ))}
                    </div>
                  </div>
                </div>
              )}
            </div>

            {/* Molecular Descriptors */}
            <div>
              <button
                onClick={() => setExpandedSection(expandedSection === 'descriptors' ? null : 'descriptors')}
                className="w-full flex justify-between items-center border-b border-journal-border pb-4 mb-4 hover:text-journal-accent transition-colors"
              >
                <div className="flex items-center gap-2">
                  <h3 className="font-bold text-xs uppercase tracking-widest text-journal-dim">2. Molecular Fingerprint (15 Descriptors)</h3>
                  <span className="text-xs text-journal-accent bg-amber-50 px-2 py-1 rounded-full">Crucial for prediction</span>
                </div>
                {expandedSection === 'descriptors' ? <ChevronDown size={20} /> : <ChevronRight size={20} />}
              </button>

              {expandedSection === 'descriptors' && (
                <div className="bg-slate-50/50 p-6 rounded-2xl border border-slate-100 animate-fade-in">
                  <div className="grid grid-cols-2 md:grid-cols-3 gap-x-6 gap-y-5">
                    {/* Primary Descriptors */}
                    <div className="col-span-2 md:col-span-3 mb-2">
                      <h4 className="text-xs font-bold text-journal-accent uppercase tracking-wider mb-3">Primary Properties</h4>
                    </div>
                    <div>
                      <label className="block text-sm font-bold text-journal-text mb-2">Molecular Weight (g/mol)</label>
                      <input type="number" step="0.1" value={formData.MolecularWeight || formData.mw || 0} onChange={(e) => updateField('MolecularWeight', parseFloat(e.target.value))} className="w-full px-3 py-2 bg-white border border-journal-border rounded-lg focus:border-journal-accent outline-none transition-all font-mono text-sm" />
                      <p className="text-[10px] text-journal-dim mt-1 italic">Size of the molecule</p>
                    </div>
                    <div>
                      <label className="block text-sm font-bold text-journal-text mb-2">LogP (Lipophilicity)</label>
                      <input type="number" step="0.1" value={formData.LogP || formData.logP || 0} onChange={(e) => updateField('LogP', parseFloat(e.target.value))} className="w-full px-3 py-2 bg-white border border-journal-border rounded-lg focus:border-journal-accent outline-none transition-all font-mono text-sm" />
                      <p className="text-[10px] text-journal-dim mt-1 italic">Fat solubility</p>
                    </div>
                    <div>
                      <label className="block text-sm font-bold text-journal-text mb-2">TPSA (Å²)</label>
                      <input type="number" step="0.1" value={formData.TPSA || 0} onChange={(e) => updateField('TPSA', parseFloat(e.target.value))} className="w-full px-3 py-2 bg-white border border-journal-border rounded-lg focus:border-journal-accent outline-none transition-all font-mono text-sm" />
                      <p className="text-[10px] text-journal-dim mt-1 italic">Polar surface area</p>
                    </div>

                    {/* Hydrogen Bonding */}
                    <div className="col-span-2 md:col-span-3 mt-4 mb-2 pt-4 border-t border-slate-200">
                      <h4 className="text-xs font-bold text-journal-accent uppercase tracking-wider mb-3">Hydrogen Bonding</h4>
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">H-Bond Donors</label>
                      <input type="number" value={formData.NumHDonors || 0} onChange={(e) => updateField('NumHDonors', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">H-Bond Acceptors</label>
                      <input type="number" value={formData.NumHAcceptors || 0} onChange={(e) => updateField('NumHAcceptors', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Rotatable Bonds</label>
                      <input type="number" value={formData.NumRotatableBonds || 0} onChange={(e) => updateField('NumRotatableBonds', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>

                    {/* Ring Systems */}
                    <div className="col-span-2 md:col-span-3 mt-4 mb-2 pt-4 border-t border-slate-200">
                      <h4 className="text-xs font-bold text-journal-accent uppercase tracking-wider mb-3">Ring Systems</h4>
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Aromatic Rings</label>
                      <input type="number" value={formData.NumAromaticRings || 0} onChange={(e) => updateField('NumAromaticRings', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Saturated Rings</label>
                      <input type="number" value={formData.NumSaturatedRings || 0} onChange={(e) => updateField('NumSaturatedRings', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Aliphatic Rings</label>
                      <input type="number" value={formData.NumAliphaticRings || 0} onChange={(e) => updateField('NumAliphaticRings', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>

                    {/* Atom Counts */}
                    <div className="col-span-2 md:col-span-3 mt-4 mb-2 pt-4 border-t border-slate-200">
                      <h4 className="text-xs font-bold text-journal-accent uppercase tracking-wider mb-3">Atom Composition</h4>
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Heteroatoms</label>
                      <input type="number" value={formData.NumHeteroatoms || 0} onChange={(e) => updateField('NumHeteroatoms', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Heavy Atoms</label>
                      <input type="number" value={formData.HeavyAtomCount || 0} onChange={(e) => updateField('HeavyAtomCount', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Fraction CSP3</label>
                      <input type="number" step="0.01" value={formData.FractionCSP3 || 0} onChange={(e) => updateField('FractionCSP3', parseFloat(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>

                    {/* Advanced Properties */}
                    <div className="col-span-2 md:col-span-3 mt-4 mb-2 pt-4 border-t border-slate-200">
                      <h4 className="text-xs font-bold text-journal-accent uppercase tracking-wider mb-3">Advanced Properties</h4>
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Molar Refractivity</label>
                      <input type="number" step="0.1" value={formData.MolarRefractivity || 0} onChange={(e) => updateField('MolarRefractivity', parseFloat(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Bertz Complexity</label>
                      <input type="number" step="0.1" value={formData.BertzCT || 0} onChange={(e) => updateField('BertzCT', parseFloat(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                    <div>
                      <label className="block text-xs font-medium text-journal-text mb-1.5">Total Rings</label>
                      <input type="number" value={formData.NumRings || 0} onChange={(e) => updateField('NumRings', parseInt(e.target.value))} className="w-full px-3 py-1.5 bg-white border border-journal-border rounded-lg outline-none text-sm font-mono" />
                    </div>
                  </div>
                </div>
              )}
            </div>

            <div className="pt-4">
              <Button onClick={handleAnalyze} disabled={loading} className="w-full text-lg">
                {loading ? <><Activity className="animate-spin" /> Analyzing...</> : 'Run Prediction Model'}
              </Button>
            </div>
          </Card>
        </div>

        {/* Output */}
        <div className="lg:col-span-5 h-full">
           <Card className={`h-full min-h-[600px] flex flex-col justify-center items-center text-center transition-all duration-gentle ${!result ? 'bg-journal-surface' : ((result.toxicity === 'Safe' || result.prediction === 0) ? 'bg-[#F0FDF4] border-[#86EFAC]' : 'bg-[#FFF1F2] border-[#FECACA]')}`}>
              {!result ? (
                <>
                   <div className="w-28 h-28 rounded-full bg-journal-bg/50 border-2 border-journal-border flex items-center justify-center mb-8">
                      <BeeIcon className="w-14 h-14 text-journal-dim opacity-15" />
                   </div>
                   <h3 className="font-display text-3xl font-bold text-journal-text mb-3">Awaiting Input</h3>
                   <p className="text-journal-dim font-serif max-w-md leading-relaxed">The bees are waiting. Enter chemical properties to assess potential toxicity risks using our gradient boosting classifier.</p>
                </>
              ) : (
                <div className="w-full h-full flex flex-col items-start text-left animate-fade-in p-2">
                   <div className="w-full flex justify-center mb-8">
                      <div className={`flex flex-col items-center ${(result.toxicity === 'Safe' || result.prediction === 0) ? 'text-journal-green' : 'text-journal-red'}`}>
                        {(result.toxicity === 'Safe' || result.prediction === 0) ? <Leaf size={64} strokeWidth={1.5} /> : <AlertTriangle size={64} strokeWidth={1.5} />}
                        <h2 className="text-5xl font-display font-bold mt-4">{(result.toxicity === 'Safe' || result.prediction === 0) ? 'Likely Safe' : 'Toxic'}</h2>
                        <p className="font-serif italic text-journal-dim mt-2">
                          {result.label_text || `Molecular profile suggests ${(result.toxicity === 'Safe' || result.prediction === 0) ? 'low' : 'high'} risk to honey bees.`}
                        </p>
                      </div>
                   </div>

                   {result.probability_toxic !== undefined && result.probability_non_toxic !== undefined ? (
                     <div className="w-full space-y-4 mb-6">
                       <div className="bg-white/70 rounded-xl p-4 border border-black/5 shadow-sm">
                         <div className="flex justify-between items-center mb-2">
                           <span className="text-xs font-bold uppercase tracking-widest text-red-600">Toxic Probability</span>
                           <span className="font-mono text-xl font-bold text-red-700">{(result.probability_toxic * 100).toFixed(1)}%</span>
                         </div>
                         <div className="w-full h-2 bg-gray-200/80 rounded-full overflow-hidden">
                           <div className="h-full bg-gradient-to-r from-red-500 to-rose-600 transition-all duration-gentle" style={{ width: `${result.probability_toxic * 100}%` }}></div>
                         </div>
                       </div>
                       <div className="bg-white/70 rounded-xl p-4 border border-black/5 shadow-sm">
                         <div className="flex justify-between items-center mb-2">
                           <span className="text-xs font-bold uppercase tracking-widest text-green-600">Non-Toxic Probability</span>
                           <span className="font-mono text-xl font-bold text-green-700">{(result.probability_non_toxic * 100).toFixed(1)}%</span>
                         </div>
                         <div className="w-full h-2 bg-gray-200/80 rounded-full overflow-hidden">
                           <div className="h-full bg-gradient-to-r from-green-500 to-emerald-600 transition-all duration-gentle" style={{ width: `${result.probability_non_toxic * 100}%` }}></div>
                         </div>
                       </div>
                     </div>
                   ) : (
                     <div className="w-full bg-white/70 rounded-xl p-6 mb-6 border border-black/5 shadow-sm">
                        <div className="flex justify-between items-end mb-3">
                           <span className="text-xs font-bold uppercase tracking-widest text-journal-dim">Model Confidence</span>
                           <span className="font-mono text-3xl font-bold">{result.confidence}%</span>
                        </div>
                        <div className="w-full h-3 bg-gray-200/80 rounded-full overflow-hidden shadow-inner">
                          <div
                            className={`h-full transition-all duration-gentle ${(result.toxicity === 'Safe' || result.prediction === 0) ? 'bg-gradient-to-r from-journal-green to-emerald-500' : 'bg-gradient-to-r from-journal-red to-rose-600'}`}
                            style={{ width: `${result.confidence}%` }}
                          ></div>
                        </div>
                        <div className="flex justify-between mt-2 text-[10px] text-gray-400 font-mono uppercase tracking-wider">
                          <span>Uncertain</span>
                          <span>Certain</span>
                        </div>
                     </div>
                   )}

                   {(result.explanation || result.recommendation) && (
                     <div className="space-y-6 w-full">
                       {result.explanation && (
                         <div>
                           <div className="flex items-center gap-3 mb-3 text-journal-accent">
                              <Microscope size={22} />
                              <h4 className="font-bold font-display text-xl">Analysis</h4>
                           </div>
                           <p className="text-journal-text leading-relaxed pl-9">{result.explanation}</p>
                         </div>
                       )}

                       {result.recommendation && (
                         <div>
                           <div className="flex items-center gap-3 mb-3 text-journal-accent">
                              <BookOpen size={22} />
                              <h4 className="font-bold font-display text-xl">Recommendation</h4>
                           </div>
                           <p className="text-journal-text leading-relaxed pl-9">{result.recommendation}</p>
                         </div>
                       )}
                     </div>
                   )}

                   <div className="mt-auto w-full pt-8 text-right">
                      <span className="font-mono text-xs text-journal-dim/50 uppercase">
                        Model: XGBoost v1.4 | {result.timestamp ? new Date(result.timestamp).toLocaleTimeString() : new Date().toLocaleTimeString()}
                      </span>
                   </div>
                </div>
              )}
           </Card>
        </div>
      </div>
    </div>
  );
};

// Scenario data from original app
const scenarioData: Scenario[] = [
  {
    id: 'almond-grower',
    title: "The Almond Grower's Decision",
    icon: '🌸',
    userType: 'Maria Rodriguez, California Almond Farmer',
    context: "Maria manages 800 acres of almonds with 1,600 rented bee colonies ($290K in rental fees). A mid-bloom aphid outbreak threatens her $8.8M crop, and she needs to spray within 48 hours while bees are actively foraging.",
    stakes: "$340K at risk (bee rentals + penalty clause)",
    compounds: [
      {
        name: 'Imidacloprid',
        description: 'Neonicotinoid insecticide - highly effective on aphids, systemic action',
        expectedResult: 'TOXIC',
        inputs: {
          source: 'PPDB',
          year: 1992,
          toxicity_type: 'Oral',
          insecticide: 1,
          herbicide: 0,
          fungicide: 0,
          other_agrochemical: 0,
          MolecularWeight: 255.66,
          LogP: 0.57,
          NumHDonors: 1,
          NumHAcceptors: 4,
          NumRotatableBonds: 2,
          NumAromaticRings: 1,
          TPSA: 86.26,
          NumHeteroatoms: 6,
          NumSaturatedRings: 1,
          NumAliphaticRings: 1,
          NumRings: 2,
          FractionCSP3: 0.25,
          MolarRefractivity: 64.2,
          BertzCT: 420.5,
          HeavyAtomCount: 17
        }
      },
      {
        name: 'Azadirachtin (Neem)',
        description: 'Botanical insecticide from neem tree - slower acting, organic-approved',
        expectedResult: 'NON-TOXIC',
        inputs: {
          source: 'BPDB',
          year: 2010,
          toxicity_type: 'Contact',
          insecticide: 1,
          herbicide: 0,
          fungicide: 0,
          other_agrochemical: 0,
          MolecularWeight: 720.71,
          LogP: 1.09,
          NumHDonors: 3,
          NumHAcceptors: 16,
          NumRotatableBonds: 7,
          NumAromaticRings: 0,
          TPSA: 215.34,
          NumHeteroatoms: 16,
          NumSaturatedRings: 6,
          NumAliphaticRings: 6,
          NumRings: 6,
          FractionCSP3: 0.76,
          MolarRefractivity: 168.5,
          BertzCT: 1850.2,
          HeavyAtomCount: 51
        }
      }
    ]
  },
  {
    id: 'rd-screening',
    title: "The R&D Startup's $3M Decision",
    icon: '🔬',
    userType: 'NexGen AgroSciences, Fungicide Development Team',
    context: "The team has synthesized 40 novel fungicide candidates. Before committing $62K per compound for EPA bee toxicity testing (OECD 213/214/245), they need to prioritize which candidates to advance.",
    stakes: "$2.5M potential savings + 6-month timeline acceleration",
    compounds: [
      {
        name: 'Candidate A (Triazole-based)',
        description: 'High LogP (4.2), multiple aromatic rings - structural alerts present',
        expectedResult: 'TOXIC',
        inputs: {
          source: 'PPDB',
          year: 2024,
          toxicity_type: 'Contact',
          insecticide: 0,
          herbicide: 0,
          fungicide: 1,
          other_agrochemical: 0,
          MolecularWeight: 385.5,
          LogP: 4.2,
          NumHDonors: 1,
          NumHAcceptors: 5,
          NumRotatableBonds: 6,
          NumAromaticRings: 3,
          TPSA: 55.8,
          NumHeteroatoms: 6,
          NumSaturatedRings: 0,
          NumAliphaticRings: 0,
          NumRings: 3,
          FractionCSP3: 0.18,
          MolarRefractivity: 105.3,
          BertzCT: 920.4,
          HeavyAtomCount: 27
        }
      },
      {
        name: 'Candidate B (Pyrimidine-based)',
        description: 'Moderate LogP (2.8), optimized structure - favorable safety profile',
        expectedResult: 'NON-TOXIC',
        inputs: {
          source: 'PPDB',
          year: 2024,
          toxicity_type: 'Contact',
          insecticide: 0,
          herbicide: 0,
          fungicide: 1,
          other_agrochemical: 0,
          MolecularWeight: 298.3,
          LogP: 2.8,
          NumHDonors: 1,
          NumHAcceptors: 4,
          NumRotatableBonds: 4,
          NumAromaticRings: 1,
          TPSA: 68.2,
          NumHeteroatoms: 5,
          NumSaturatedRings: 1,
          NumAliphaticRings: 1,
          NumRings: 2,
          FractionCSP3: 0.42,
          MolarRefractivity: 78.6,
          BertzCT: 580.2,
          HeavyAtomCount: 21
        }
      }
    ]
  },
  {
    id: 'beekeeper-investigation',
    title: "The Beekeeper's Investigation",
    icon: '🔍',
    userType: 'James Patterson, Commercial Beekeeper (Iowa)',
    context: "James discovered 48 of his 200 colonies collapsed with dead bees at hive entrances. He suspects pesticide exposure from neighboring corn/soybean farms but lab results take 2-3 weeks. He needs to decide whether to move remaining colonies immediately.",
    stakes: "$52K at risk (lost colonies + prevented losses)",
    compounds: [
      {
        name: 'Clothianidin',
        description: 'Neonicotinoid - common corn seed treatment, systemic',
        expectedResult: 'TOXIC',
        inputs: {
          source: 'PPDB',
          year: 2001,
          toxicity_type: 'Oral',
          insecticide: 1,
          herbicide: 0,
          fungicide: 0,
          other_agrochemical: 0,
          MolecularWeight: 249.68,
          LogP: 0.91,
          NumHDonors: 2,
          NumHAcceptors: 4,
          NumRotatableBonds: 3,
          NumAromaticRings: 1,
          TPSA: 101.15,
          NumHeteroatoms: 7,
          NumSaturatedRings: 0,
          NumAliphaticRings: 0,
          NumRings: 1,
          FractionCSP3: 0.17,
          MolarRefractivity: 62.8,
          BertzCT: 385.2,
          HeavyAtomCount: 16
        }
      },
      {
        name: 'Glyphosate',
        description: 'Herbicide - targets plant enzyme (EPSPS) that bees lack',
        expectedResult: 'NON-TOXIC',
        inputs: {
          source: 'PPDB',
          year: 1999,
          toxicity_type: 'Contact',
          insecticide: 0,
          herbicide: 1,
          fungicide: 0,
          other_agrochemical: 0,
          MolecularWeight: 169.07,
          LogP: -3.2,
          NumHDonors: 4,
          NumHAcceptors: 6,
          NumRotatableBonds: 4,
          NumAromaticRings: 0,
          TPSA: 108.21,
          NumHeteroatoms: 6,
          NumSaturatedRings: 0,
          NumAliphaticRings: 0,
          NumRings: 0,
          FractionCSP3: 0.67,
          MolarRefractivity: 35.2,
          BertzCT: 125.8,
          HeavyAtomCount: 12
        }
      },
      {
        name: 'Atrazine',
        description: 'Herbicide - common corn herbicide, triazine class',
        expectedResult: 'NON-TOXIC',
        inputs: {
          source: 'ECOTOX',
          year: 1961,
          toxicity_type: 'Contact',
          insecticide: 0,
          herbicide: 1,
          fungicide: 0,
          other_agrochemical: 0,
          MolecularWeight: 215.68,
          LogP: 2.61,
          NumHDonors: 2,
          NumHAcceptors: 5,
          NumRotatableBonds: 4,
          NumAromaticRings: 1,
          TPSA: 62.73,
          NumHeteroatoms: 6,
          NumSaturatedRings: 0,
          NumAliphaticRings: 0,
          NumRings: 1,
          FractionCSP3: 0.5,
          MolarRefractivity: 58.9,
          BertzCT: 245.6,
          HeavyAtomCount: 14
        }
      }
    ]
  }
];

const ScenariosView = () => {
  const [activeScenario, setActiveScenario] = useState<string>(scenarioData[0].id);
  const [selectedCompound, setSelectedCompound] = useState<number>(0);
  const [results, setResults] = useState<Record<string, PredictionResult | null>>({});
  const [loading, setLoading] = useState<string | null>(null);
  const [editedInputs, setEditedInputs] = useState<Record<string, Partial<ChemicalData>>>({});

  const currentScenario = scenarioData.find(s => s.id === activeScenario)!;
  const currentCompound = currentScenario.compounds[selectedCompound];

  const getInputsKey = (scenarioId: string, compoundIdx: number) => `${scenarioId}-${compoundIdx}`;
  const currentInputs = editedInputs[getInputsKey(activeScenario, selectedCompound)] || currentCompound.inputs;

  const handleInputChange = (field: string, value: number | string) => {
    const key = getInputsKey(activeScenario, selectedCompound);
    setEditedInputs(prev => ({
      ...prev,
      [key]: {
        ...(prev[key] || currentCompound.inputs),
        [field]: value
      }
    }));
  };

  const runPrediction = async (compound: CompoundData, index: number) => {
    const key = getInputsKey(activeScenario, index);
    const inputs = editedInputs[key] || compound.inputs;
    setLoading(key);

    try {
      const result = await analyzeChemicalToxicity(inputs as ChemicalData);
      setResults(prev => ({ ...prev, [key]: result }));
    } catch (err) {
      console.error('Prediction failed:', err);
    } finally {
      setLoading(null);
    }
  };

  const runAllPredictions = async () => {
    for (let i = 0; i < currentScenario.compounds.length; i++) {
      await runPrediction(currentScenario.compounds[i], i);
    }
  };

  const getResultForCompound = (index: number) => {
    const key = getInputsKey(activeScenario, index);
    return results[key];
  };

  return (
    <div className="animate-fade-in space-y-6">
      <SectionHeader
        title="Real-World Scenarios"
        subtitle="Explore how ApisTox helps different stakeholders make data-driven decisions about pesticide safety."
      />

      {/* Scenario Selector */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        {scenarioData.map((scenario) => (
          <button
            key={scenario.id}
            onClick={() => {
              setActiveScenario(scenario.id);
              setSelectedCompound(0);
            }}
            className={`p-6 border-2 rounded-xl text-left transition-all shadow-card-subtle hover:shadow-lifted ${
              activeScenario === scenario.id
                ? 'border-journal-accent bg-[#FEF3C7] ring-2 ring-journal-accent/20'
                : 'border-journal-border bg-white hover:border-journal-accent/40'
            }`}
          >
            <div className="text-5xl mb-4">{scenario.icon}</div>
            <h3 className="font-bold font-display text-lg mb-2">{scenario.title}</h3>
            <p className="text-xs font-bold text-journal-dim uppercase tracking-wider">{scenario.userType.split(',')[0]}</p>
          </button>
        ))}
      </div>

      {/* Active Scenario Detail */}
      <Card className="p-0 overflow-hidden">
        {/* Scenario Header */}
        <div className="bg-gradient-to-r from-[#FEF3C7] to-[#FED7AA] p-8 border-b border-journal-border">
          <div className="flex items-start gap-4">
            <div className="text-5xl">{currentScenario.icon}</div>
            <div className="flex-1">
              <h3 className="text-2xl font-bold font-display text-journal-text">{currentScenario.title}</h3>
              <p className="text-journal-maroon font-medium mt-1">{currentScenario.userType}</p>
              <p className="text-journal-text mt-3 leading-relaxed">{currentScenario.context}</p>
              <div className="mt-4 inline-flex items-center bg-journal-accent/15 text-journal-accent px-4 py-2 rounded-full text-sm font-bold">
                💰 {currentScenario.stakes}
              </div>
            </div>
          </div>
        </div>

        {/* Compound Tabs */}
        <div className="border-b border-gray-100 px-6 pt-4 bg-white">
          <div className="flex gap-2 overflow-x-auto pb-4">
            {currentScenario.compounds.map((compound, idx) => {
              const result = getResultForCompound(idx);
              return (
                <button
                  key={idx}
                  onClick={() => setSelectedCompound(idx)}
                  className={`px-4 py-2 rounded-lg font-medium text-sm whitespace-nowrap transition-all flex items-center gap-2 ${
                    selectedCompound === idx
                      ? 'bg-journal-accent text-white shadow-card-subtle'
                      : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                  }`}
                >
                  {compound.name}
                  {result && (
                    <span className={`w-2 h-2 rounded-full ${
                      (result.prediction === 1 || result.toxicity === 'Toxic') ? 'bg-red-400' : 'bg-green-400'
                    }`} />
                  )}
                </button>
              );
            })}
            <button
              onClick={runAllPredictions}
              className="px-4 py-2 rounded-lg font-medium text-sm whitespace-nowrap bg-journal-maroon text-white hover:bg-[#5C1D0A] transition-all flex items-center gap-2 ml-auto"
            >
              Run All <Play size={16} fill="currentColor" />
            </button>
          </div>
        </div>

        {/* Compound Details + Results */}
        <div className="p-8">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            {/* Left: Compound Info + Editable Inputs */}
            <div className="space-y-4">
              <div className="bg-gray-50 rounded-xl p-6 border border-gray-200">
                <h4 className="font-bold font-display text-xl text-journal-text">{currentCompound.name}</h4>
                <p className="text-sm text-journal-dim mt-2">{currentCompound.description}</p>
                <div className="mt-3 flex items-center gap-2">
                  <span className="text-xs text-gray-500">Expected:</span>
                  <span className={`text-xs font-bold px-3 py-1 rounded-full ${
                    currentCompound.expectedResult === 'TOXIC'
                      ? 'bg-red-100 text-red-700'
                      : 'bg-green-100 text-green-700'
                  }`}>
                    {currentCompound.expectedResult}
                  </span>
                </div>
              </div>

              {/* Key Editable Fields */}
              <div className="space-y-3 bg-slate-50 p-4 rounded-xl border border-slate-200">
                <h5 className="font-medium text-journal-text text-sm mb-3">Key Parameters (editable)</h5>

                <div className="grid grid-cols-2 gap-3">
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">Year</label>
                    <input
                      type="number"
                      value={currentInputs.year || 2020}
                      onChange={(e) => handleInputChange('year', parseInt(e.target.value))}
                      className="w-full px-3 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">Source</label>
                    <select
                      value={currentInputs.source || 'PPDB'}
                      onChange={(e) => handleInputChange('source', e.target.value)}
                      className="w-full px-3 py-2 text-sm border border-gray-200 rounded-lg bg-white focus:border-journal-accent outline-none"
                    >
                      <option value="PPDB">PPDB</option>
                      <option value="ECOTOX">ECOTOX</option>
                      <option value="BPDB">BPDB</option>
                    </select>
                  </div>
                </div>

                <div className="grid grid-cols-3 gap-3">
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">Mol. Wt</label>
                    <input
                      type="number"
                      step="0.1"
                      value={currentInputs.MolecularWeight || 0}
                      onChange={(e) => handleInputChange('MolecularWeight', parseFloat(e.target.value))}
                      className="w-full px-2 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">LogP</label>
                    <input
                      type="number"
                      step="0.1"
                      value={currentInputs.LogP || 0}
                      onChange={(e) => handleInputChange('LogP', parseFloat(e.target.value))}
                      className="w-full px-2 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">TPSA</label>
                    <input
                      type="number"
                      step="0.1"
                      value={currentInputs.TPSA || 0}
                      onChange={(e) => handleInputChange('TPSA', parseFloat(e.target.value))}
                      className="w-full px-2 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                </div>

                <div className="grid grid-cols-3 gap-3">
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">Aromatic Rings</label>
                    <input
                      type="number"
                      value={currentInputs.NumAromaticRings || 0}
                      onChange={(e) => handleInputChange('NumAromaticRings', parseInt(e.target.value))}
                      className="w-full px-2 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">H-Donors</label>
                    <input
                      type="number"
                      value={currentInputs.NumHDonors || 0}
                      onChange={(e) => handleInputChange('NumHDonors', parseInt(e.target.value))}
                      className="w-full px-2 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                  <div>
                    <label className="text-xs text-gray-600 font-medium block mb-1">H-Acceptors</label>
                    <input
                      type="number"
                      value={currentInputs.NumHAcceptors || 0}
                      onChange={(e) => handleInputChange('NumHAcceptors', parseInt(e.target.value))}
                      className="w-full px-2 py-2 text-sm border border-gray-200 rounded-lg bg-white font-mono focus:border-journal-accent outline-none"
                    />
                  </div>
                </div>
              </div>

              <Button onClick={() => runPrediction(currentCompound, selectedCompound)} disabled={loading !== null} className="w-full">
                {loading === getInputsKey(activeScenario, selectedCompound) ? (
                  <><Activity className="animate-spin" /> Analyzing...</>
                ) : (
                  `Run Analysis for ${currentCompound.name}`
                )}
              </Button>
            </div>

            {/* Right: Results */}
            <div>
              {(() => {
                const result = getResultForCompound(selectedCompound);
                if (!result) {
                  return (
                    <div className="h-full flex items-center justify-center bg-gray-50 rounded-xl border-2 border-dashed border-gray-200 p-8">
                      <div className="text-center text-gray-500">
                        <div className="text-4xl mb-3">🔮</div>
                        <p className="font-medium">Run analysis to see prediction</p>
                        <p className="text-sm mt-1">Results will appear here</p>
                      </div>
                    </div>
                  );
                }

                const isToxic = result.prediction === 1 || result.toxicity === 'Toxic';
                const probToxic = result.probability_toxic ?? (isToxic ? 0.85 : 0.15);
                const probNonToxic = result.probability_non_toxic ?? (isToxic ? 0.15 : 0.85);

                return (
                  <div className={`rounded-xl p-6 ${isToxic ? 'bg-red-50 border-2 border-red-200' : 'bg-green-50 border-2 border-green-200'}`}>
                    <div className="text-center mb-6">
                      <div className="text-6xl mb-3">{isToxic ? '⚠️' : '✅'}</div>
                      <div className={`text-3xl font-bold font-display ${isToxic ? 'text-red-700' : 'text-green-700'}`}>
                        {isToxic ? 'TOXIC' : 'NON-TOXIC'}
                      </div>
                      <div className="text-sm text-gray-600 mt-1">to honey bees (Apis mellifera)</div>
                    </div>

                    <div className="space-y-3">
                      <div className="bg-white rounded-lg p-3">
                        <div className="flex justify-between text-sm mb-1">
                          <span className="text-gray-600">Toxic Probability</span>
                          <span className="font-bold text-red-600">{(probToxic * 100).toFixed(1)}%</span>
                        </div>
                        <div className="h-2 bg-gray-200 rounded-full overflow-hidden">
                          <div
                            className="h-full bg-red-500 rounded-full transition-all"
                            style={{ width: `${probToxic * 100}%` }}
                          />
                        </div>
                      </div>

                      <div className="bg-white rounded-lg p-3">
                        <div className="flex justify-between text-sm mb-1">
                          <span className="text-gray-600">Non-Toxic Probability</span>
                          <span className="font-bold text-green-600">{(probNonToxic * 100).toFixed(1)}%</span>
                        </div>
                        <div className="h-2 bg-gray-200 rounded-full overflow-hidden">
                          <div
                            className="h-full bg-green-500 rounded-full transition-all"
                            style={{ width: `${probNonToxic * 100}%` }}
                          />
                        </div>
                      </div>

                      <div className="bg-white rounded-lg p-3">
                        <div className="flex justify-between text-sm">
                          <span className="text-gray-600">Model Confidence</span>
                          <span className="font-bold text-gray-800">{(result.confidence || 75).toFixed(1)}%</span>
                        </div>
                      </div>
                    </div>

                    {/* Match Check */}
                    <div className={`mt-4 p-3 rounded-lg ${
                      (isToxic && currentCompound.expectedResult === 'TOXIC') ||
                      (!isToxic && currentCompound.expectedResult === 'NON-TOXIC')
                        ? 'bg-blue-100 text-blue-800'
                        : 'bg-yellow-100 text-yellow-800'
                    }`}>
                      <div className="text-sm font-medium">
                        {(isToxic && currentCompound.expectedResult === 'TOXIC') ||
                         (!isToxic && currentCompound.expectedResult === 'NON-TOXIC')
                          ? '✓ Prediction matches expected result'
                          : '⚡ Prediction differs from expected - try adjusting parameters'
                        }
                      </div>
                    </div>
                  </div>
                );
              })()}
            </div>
          </div>
        </div>

        {/* Scenario Insight */}
        <div className="bg-slate-50 p-6 border-t border-slate-200">
          <h4 className="font-bold text-slate-800 mb-2 flex items-center gap-2"><Info size={18} /> Key Insight</h4>
          {activeScenario === 'almond-grower' && (
            <p className="text-slate-600 text-sm leading-relaxed">
              <strong>Neonicotinoids are systemic insecticides</strong> that persist in plant tissue and nectar, making them highly toxic to bees even days after application.
              Botanical alternatives like azadirachtin break down faster and have lower acute toxicity, though they may require more frequent application.
              <strong> The "insecticide" flag is the #1 predictor of bee toxicity in our SHAP analysis.</strong>
            </p>
          )}
          {activeScenario === 'rd-screening' && (
            <p className="text-slate-600 text-sm leading-relaxed">
              <strong>High LogP (&gt;3.5) combined with multiple aromatic rings</strong> correlates with increased bee toxicity due to bioaccumulation potential.
              Our model identifies these structural alerts early, allowing R&D teams to optimize molecular design before expensive regulatory testing.
              <strong> Screening 40 compounds takes minutes vs. 18+ months and $2.5M for traditional OECD testing.</strong>
            </p>
          )}
          {activeScenario === 'beekeeper-investigation' && (
            <p className="text-slate-600 text-sm leading-relaxed">
              <strong>Herbicides target plant-specific enzymes (like EPSPS for glyphosate)</strong> that bees don't possess, making them generally non-toxic to pollinators.
              Neonicotinoid seed treatments, however, are a leading cause of bee kills - the active ingredient persists in soil and is taken up by subsequent crops.
              <strong> ApisTox helps prioritize investigation targets before lab results arrive.</strong>
            </p>
          )}
        </div>
      </Card>
    </div>
  );
};

type ScienceTabId = 'data' | 'model' | 'science';

const ScienceView = () => {
  const [activeSubTab, setActiveSubTab] = useState<ScienceTabId>('data');

  const tabs = [
    { id: 'data' as ScienceTabId, label: 'The Data', icon: '📊' },
    { id: 'model' as ScienceTabId, label: 'The Model', icon: '⚙️' },
    { id: 'science' as ScienceTabId, label: 'The Science', icon: '🔬' },
  ];

  return (
    <div className="animate-fade-in space-y-6">
      <SectionHeader
        title="Scientific Basis"
        subtitle="Understanding the data, model architecture, and interpretability behind ApisTox predictions."
      />

      {/* Sub-Tab Navigation */}
      <div className="bg-white rounded-2xl shadow-card-subtle border border-journal-border p-2">
        <div className="flex space-x-2">
          {tabs.map((tab) => (
            <button
              key={tab.id}
              onClick={() => setActiveSubTab(tab.id)}
              className={`flex-1 py-3 px-4 rounded-xl font-medium transition-all duration-quick ${
                activeSubTab === tab.id
                  ? 'bg-gradient-to-r from-journal-accent to-amber-600 text-white shadow-card-subtle'
                  : 'text-journal-dim hover:bg-amber-50'
              }`}
            >
              <span className="mr-2">{tab.icon}</span>
              {tab.label}
            </button>
          ))}
        </div>
      </div>

      {/* Sub-Tab Content */}
      <Card className="min-h-[600px]">
        {activeSubTab === 'data' && <DataTab />}
        {activeSubTab === 'model' && <ModelTab />}
        {activeSubTab === 'science' && <ScienceTab />}
      </Card>
    </div>
  );
};

// Data Tab Component
const DataTab = () => {
  return (
    <div className="space-y-10">
      <div className="text-center max-w-2xl mx-auto">
        <h2 className="text-3xl font-bold font-display text-journal-text mb-4">Understanding the Data</h2>
        <p className="text-journal-dim font-serif leading-relaxed">
          Our model learns from real-world toxicity data collected from three major scientific databases. Here's what powers our predictions.
        </p>
      </div>

      {/* What We're Predicting */}
      <section>
        <div className="flex items-center gap-3 mb-6">
          <span className="w-10 h-10 rounded-full bg-journal-accent/20 text-journal-accent font-bold flex items-center justify-center font-mono">1</span>
          <h3 className="text-2xl font-bold font-display">What We're Predicting</h3>
        </div>
        <div className="grid md:grid-cols-2 gap-6 bg-[#FAFAF8] p-8 rounded-2xl border-2 border-gray-100">
          <div className="bg-white p-8 rounded-xl text-center shadow-card-subtle border-2 border-emerald-100">
            <div className="inline-flex p-4 rounded-full bg-emerald-100 text-emerald-600 mb-4">
              <CheckCircle2 size={36} />
            </div>
            <h4 className="text-2xl font-bold font-display text-emerald-800 mb-2">Non-Toxic</h4>
            <p className="font-mono text-emerald-600 mb-2 text-lg">LD50 &gt; 11 μg/bee</p>
            <p className="text-sm text-gray-500 font-serif italic">Safe for pollinators</p>
          </div>
          <div className="bg-white p-8 rounded-xl text-center shadow-card-subtle border-2 border-rose-100">
            <div className="inline-flex p-4 rounded-full bg-rose-100 text-rose-600 mb-4">
              <AlertTriangle size={36} />
            </div>
            <h4 className="text-2xl font-bold font-display text-rose-800 mb-2">Toxic</h4>
            <p className="font-mono text-rose-600 mb-2 text-lg">LD50 ≤ 11 μg/bee</p>
            <p className="text-sm text-gray-500 font-serif italic">Harmful to honey bees</p>
          </div>
          <div className="md:col-span-2 text-center pt-4">
            <p className="text-sm font-bold text-journal-dim">LD50 = Lethal Dose that kills 50% of test population (EPA standard threshold)</p>
          </div>
        </div>
      </section>

      {/* Data Sources */}
      <section>
        <div className="flex items-center gap-3 mb-6">
          <span className="w-10 h-10 rounded-full bg-journal-accent/20 text-journal-accent font-bold flex items-center justify-center font-mono">2</span>
          <h3 className="text-2xl font-bold font-display">Data Sources</h3>
        </div>
        <div className="grid md:grid-cols-3 gap-6">
          <div className="p-6 bg-blue-50 rounded-xl border-2 border-blue-100 shadow-card-subtle">
            <h4 className="font-bold font-display text-blue-900 mb-2 text-xl">ECOTOX</h4>
            <p className="text-xs font-bold uppercase text-blue-600 mb-3 tracking-wider">EPA Ecotoxicology Database</p>
            <p className="text-sm text-blue-800 leading-relaxed">
              U.S. Environmental Protection Agency's curated database of chemical toxicity studies on aquatic and terrestrial species.
            </p>
          </div>
          <div className="p-6 bg-green-50 rounded-xl border-2 border-green-100 shadow-card-subtle">
            <h4 className="font-bold font-display text-green-900 mb-2 text-xl">PPDB</h4>
            <p className="text-xs font-bold uppercase text-green-600 mb-3 tracking-wider">Pesticide Properties Database</p>
            <p className="text-sm text-green-800 leading-relaxed">
              University of Hertfordshire's comprehensive database covering physicochemical and ecotoxicological properties.
            </p>
          </div>
          <div className="p-6 bg-purple-50 rounded-xl border-2 border-purple-100 shadow-card-subtle">
            <h4 className="font-bold font-display text-purple-900 mb-2 text-xl">BPDB</h4>
            <p className="text-xs font-bold uppercase text-purple-600 mb-3 tracking-wider">Bio-Pesticides Database</p>
            <p className="text-sm text-purple-800 leading-relaxed">
              Specialized database for biological pesticides and their environmental impact profiles.
            </p>
          </div>
        </div>
      </section>

      {/* Dataset Statistics */}
      <section>
        <div className="flex items-center gap-3 mb-6">
          <span className="w-10 h-10 rounded-full bg-journal-accent/20 text-journal-accent font-bold flex items-center justify-center font-mono">3</span>
          <h3 className="text-2xl font-bold font-display">Dataset at a Glance</h3>
        </div>
        <div className="grid md:grid-cols-2 gap-6">
          <div className="space-y-4">
            <div className="bg-gray-50 rounded-xl p-5 border border-gray-200">
              <div className="flex justify-between items-center">
                <span className="text-gray-600 font-medium">Total Compounds</span>
                <span className="text-4xl font-bold font-display text-journal-text">1,035</span>
              </div>
            </div>
            <div className="bg-green-50 rounded-xl p-5 border border-green-200">
              <div className="flex justify-between items-center mb-2">
                <span className="text-green-700 font-medium">Non-Toxic (Class 0)</span>
                <span className="text-2xl font-bold text-green-800">739 <span className="text-sm font-normal">(71.4%)</span></span>
              </div>
              <div className="bg-green-200 rounded-full h-2">
                <div className="bg-green-500 h-2 rounded-full" style={{ width: '71.4%' }}></div>
              </div>
            </div>
            <div className="bg-red-50 rounded-xl p-5 border border-red-200">
              <div className="flex justify-between items-center mb-2">
                <span className="text-red-700 font-medium">Toxic (Class 1)</span>
                <span className="text-2xl font-bold text-red-800">296 <span className="text-sm font-normal">(28.6%)</span></span>
              </div>
              <div className="bg-red-200 rounded-full h-2">
                <div className="bg-red-500 h-2 rounded-full" style={{ width: '28.6%' }}></div>
              </div>
            </div>
            <div className="bg-amber-50 rounded-xl p-4 border border-amber-200">
              <p className="text-sm text-amber-800">
                <strong>Class Imbalance:</strong> 2.5:1 ratio — we use SMOTE resampling to balance training data.
              </p>
            </div>
          </div>

          <div className="bg-gray-50 rounded-xl p-6 border border-gray-200 flex items-center justify-center">
            <div className="text-center text-gray-400">
              <BarChart3 size={64} className="mx-auto mb-3 opacity-30" />
              <p className="text-sm font-medium">Target Distribution Chart</p>
              <p className="text-xs mt-1">(Figure placeholder)</p>
            </div>
          </div>
        </div>
      </section>

      {/* No Missing Values */}
      <section className="bg-gradient-to-r from-emerald-50 to-teal-50 rounded-xl p-6 border border-emerald-200">
        <div className="flex items-start space-x-4">
          <div className="text-3xl">✨</div>
          <div>
            <h4 className="font-bold text-emerald-800 mb-1">Clean Data, No Missing Values</h4>
            <p className="text-sm text-emerald-700 leading-relaxed">
              The ApisTox dataset is fully curated with no missing values. Every compound has complete molecular structure (SMILES),
              toxicity measurements, and metadata — ready for machine learning without imputation.
            </p>
          </div>
        </div>
      </section>
    </div>
  );
};

// Model Tab Component
const ModelTab = () => {
  return (
    <div className="space-y-10">
      <div className="text-center max-w-2xl mx-auto">
        <h2 className="text-3xl font-bold font-display text-journal-text mb-4">How the Model Works</h2>
        <p className="text-journal-dim font-serif leading-relaxed">
          From raw chemical data to accurate predictions — here's the machine learning pipeline that powers ApisTox.
        </p>
      </div>

      {/* Pipeline Flowchart */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">ML Pipeline Overview</h3>
        <div className="bg-gradient-to-r from-slate-50 to-slate-100 rounded-xl p-6 border border-slate-200">
          <div className="flex flex-wrap justify-center items-center gap-3 text-sm">
            <div className="bg-blue-500 text-white px-4 py-2 rounded-lg font-medium shadow-card-subtle">📥 Data Ingestion</div>
            <div className="text-slate-400">→</div>
            <div className="bg-purple-500 text-white px-4 py-2 rounded-lg font-medium shadow-card-subtle">🧬 Feature Engineering</div>
            <div className="text-slate-400">→</div>
            <div className="bg-amber-500 text-white px-4 py-2 rounded-lg font-medium shadow-card-subtle">⚙️ Preprocessing</div>
            <div className="text-slate-400">→</div>
            <div className="bg-green-500 text-white px-4 py-2 rounded-lg font-medium shadow-card-subtle">🤖 Model Training</div>
            <div className="text-slate-400">→</div>
            <div className="bg-red-500 text-white px-4 py-2 rounded-lg font-medium shadow-card-subtle">📊 Evaluation</div>
          </div>
        </div>
      </section>

      {/* Feature Engineering */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">Feature Engineering: 15 Molecular Descriptors</h3>
        <p className="text-journal-dim mb-4">
          We extract chemical properties from SMILES notation using RDKit. These descriptors capture the molecular characteristics that influence toxicity.
        </p>

        <div className="grid md:grid-cols-2 gap-6">
          <div className="bg-white rounded-xl border border-gray-200 overflow-hidden shadow-card-subtle">
            <table className="w-full text-sm">
              <thead className="bg-gray-50">
                <tr>
                  <th className="px-4 py-3 text-left font-medium text-gray-700">Descriptor</th>
                  <th className="px-4 py-3 text-left font-medium text-gray-700">Meaning</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-100">
                <tr><td className="px-4 py-2 font-mono text-xs">MolecularWeight</td><td className="px-4 py-2 text-gray-600">Mass in Daltons</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">LogP</td><td className="px-4 py-2 text-gray-600">Lipophilicity (fat solubility)</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">NumHDonors</td><td className="px-4 py-2 text-gray-600">Hydrogen bond donors</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">NumHAcceptors</td><td className="px-4 py-2 text-gray-600">Hydrogen bond acceptors</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">TPSA</td><td className="px-4 py-2 text-gray-600">Polar surface area</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">NumRotatableBonds</td><td className="px-4 py-2 text-gray-600">Molecular flexibility</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">NumAromaticRings</td><td className="px-4 py-2 text-gray-600">Aromatic ring count</td></tr>
                <tr><td className="px-4 py-2 font-mono text-xs">HeavyAtomCount</td><td className="px-4 py-2 text-gray-600">Non-hydrogen atoms</td></tr>
              </tbody>
            </table>
            <div className="px-4 py-2 bg-gray-50 text-xs text-gray-500">+ 7 more descriptors</div>
          </div>

          <div className="bg-gray-50 rounded-xl p-6 border border-gray-200 flex items-center justify-center">
            <div className="text-center text-gray-400">
              <Atom size={64} className="mx-auto mb-3 opacity-30" />
              <p className="text-sm font-medium">Molecular Descriptor Distributions</p>
              <p className="text-xs mt-1">(Figure placeholder)</p>
            </div>
          </div>
        </div>
      </section>

      {/* Preprocessing Steps */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">Preprocessing Pipeline</h3>
        <div className="grid md:grid-cols-3 gap-4">
          <div className="bg-white rounded-xl p-5 border border-gray-200 shadow-card-subtle">
            <div className="text-2xl mb-2">🔤</div>
            <h4 className="font-bold text-journal-text mb-2">One-Hot Encoding</h4>
            <p className="text-sm text-journal-dim">
              Categorical features (source, toxicity_type) converted to binary columns for model compatibility.
            </p>
          </div>
          <div className="bg-white rounded-xl p-5 border border-gray-200 shadow-card-subtle">
            <div className="text-2xl mb-2">📏</div>
            <h4 className="font-bold text-journal-text mb-2">StandardScaler</h4>
            <p className="text-sm text-journal-dim">
              Numerical features normalized to zero mean and unit variance. Fit on training data only (no leakage).
            </p>
          </div>
          <div className="bg-white rounded-xl p-5 border border-gray-200 shadow-card-subtle">
            <div className="text-2xl mb-2">⚖️</div>
            <h4 className="font-bold text-journal-text mb-2">SMOTE Resampling</h4>
            <p className="text-sm text-journal-dim">
              Synthetic Minority Over-sampling balances the 2.5:1 class imbalance by generating synthetic toxic samples.
            </p>
          </div>
        </div>
      </section>

      {/* Model Comparison */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">Model Comparison</h3>
        <div className="bg-white rounded-xl border border-gray-200 overflow-hidden shadow-card-subtle">
          <table className="w-full text-sm">
            <thead className="bg-slate-800 text-white">
              <tr>
                <th className="px-4 py-3 text-left">Model</th>
                <th className="px-4 py-3 text-center">Accuracy</th>
                <th className="px-4 py-3 text-center">F1 Score</th>
                <th className="px-4 py-3 text-center">ROC-AUC</th>
                <th className="px-4 py-3 text-center">Status</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-100">
              <tr className="bg-amber-50">
                <td className="px-4 py-3 font-bold text-amber-800">XGBoost</td>
                <td className="px-4 py-3 text-center font-mono">83.6%</td>
                <td className="px-4 py-3 text-center font-mono">70.2%</td>
                <td className="px-4 py-3 text-center font-mono">85.8%</td>
                <td className="px-4 py-3 text-center"><span className="bg-amber-500 text-white text-xs px-2 py-1 rounded-full">Selected</span></td>
              </tr>
              <tr>
                <td className="px-4 py-3 text-gray-700">Random Forest</td>
                <td className="px-4 py-3 text-center font-mono text-gray-600">81.2%</td>
                <td className="px-4 py-3 text-center font-mono text-gray-600">67.5%</td>
                <td className="px-4 py-3 text-center font-mono text-gray-600">83.1%</td>
                <td className="px-4 py-3 text-center text-gray-400">—</td>
              </tr>
              <tr>
                <td className="px-4 py-3 text-gray-700">Logistic Regression</td>
                <td className="px-4 py-3 text-center font-mono text-gray-600">78.4%</td>
                <td className="px-4 py-3 text-center font-mono text-gray-600">62.1%</td>
                <td className="px-4 py-3 text-center font-mono text-gray-600">79.5%</td>
                <td className="px-4 py-3 text-center text-gray-400">—</td>
              </tr>
            </tbody>
          </table>
        </div>
        <p className="text-sm text-journal-dim mt-3">
          XGBoost selected based on highest validation F1 score. Gradient boosting excels at structured tabular data with feature interactions.
        </p>
      </section>
    </div>
  );
};

// Science Tab Component
const ScienceTab = () => {
  return (
    <div className="space-y-10">
      <div className="text-center max-w-2xl mx-auto">
        <h2 className="text-3xl font-bold font-display text-journal-text mb-4">Deep Dive: Model Interpretability</h2>
        <p className="text-journal-dim font-serif leading-relaxed">
          Understanding <em>why</em> the model makes predictions is as important as accuracy.
          We use SHAP (SHapley Additive exPlanations) to explain feature contributions.
        </p>
      </div>

      {/* SHAP Summary Plot */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">SHAP Summary Plot</h3>
        <div className="bg-white rounded-xl p-6 border border-gray-200 shadow-card-subtle">
          <div className="flex items-center justify-center h-64 bg-gray-50 rounded-lg">
            <div className="text-center text-gray-400">
              <FileText size={64} className="mx-auto mb-3 opacity-30" />
              <p className="text-sm font-medium">SHAP Beeswarm Plot</p>
              <p className="text-xs mt-1">(Figure placeholder: /figures/shap_summary.png)</p>
            </div>
          </div>
        </div>
        <div className="mt-4 bg-blue-50 rounded-xl p-5 border border-blue-200">
          <h4 className="font-bold text-blue-800 mb-2">How to Read This Plot</h4>
          <ul className="text-sm text-blue-700 space-y-1">
            <li>• Each dot represents one compound prediction</li>
            <li>• X-axis: SHAP value (positive = pushes toward toxic, negative = pushes toward non-toxic)</li>
            <li>• Color: Feature value (red = high, blue = low)</li>
            <li>• Features sorted by importance (most impactful at top)</li>
          </ul>
        </div>
      </section>

      {/* Feature Importance */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">Feature Importance Ranking</h3>
        <div className="grid md:grid-cols-2 gap-6">
          <div className="bg-white rounded-xl p-6 border border-gray-200 shadow-card-subtle">
            <div className="flex items-center justify-center h-64 bg-gray-50 rounded-lg">
              <div className="text-center text-gray-400">
                <BarChart3 size={64} className="mx-auto mb-3 opacity-30" />
                <p className="text-sm font-medium">SHAP Feature Importance</p>
                <p className="text-xs mt-1">(Figure placeholder)</p>
              </div>
            </div>
          </div>
          <div className="space-y-4">
            <h4 className="font-bold text-journal-text">Top Predictors Explained</h4>
            <div className="space-y-3">
              <div className="bg-red-50 rounded-lg p-4 border-l-4 border-red-500">
                <div className="font-bold text-red-800">1. Insecticide Flag</div>
                <p className="text-sm text-red-700 mt-1">
                  Compounds designed to kill insects are inherently more toxic to bees (also insects). Strongest predictor.
                </p>
              </div>
              <div className="bg-orange-50 rounded-lg p-4 border-l-4 border-orange-500">
                <div className="font-bold text-orange-800">2. Herbicide Flag</div>
                <p className="text-sm text-orange-700 mt-1">
                  Herbicides target plants, not insects — generally safer for pollinators.
                </p>
              </div>
              <div className="bg-amber-50 rounded-lg p-4 border-l-4 border-amber-500">
                <div className="font-bold text-amber-800">3. Year of Registration</div>
                <p className="text-sm text-amber-700 mt-1">
                  Newer pesticides (post-2000) tend to be less toxic due to stricter environmental regulations.
                </p>
              </div>
              <div className="bg-yellow-50 rounded-lg p-4 border-l-4 border-yellow-500">
                <div className="font-bold text-yellow-800">4. LogP (Lipophilicity)</div>
                <p className="text-sm text-yellow-700 mt-1">
                  Higher fat solubility correlates with bioaccumulation and increased toxicity.
                </p>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Limitations & Ethics */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">Limitations & Ethical Use</h3>
        <div className="grid md:grid-cols-2 gap-4">
          <div className="bg-slate-50 rounded-xl p-5 border border-slate-200">
            <h4 className="font-bold text-slate-800 mb-3">⚠️ Model Limitations</h4>
            <ul className="text-sm text-slate-600 space-y-2">
              <li>• Predictions are probabilistic, not definitive</li>
              <li>• Training data may not cover all chemical classes</li>
              <li>• Does not account for synergistic effects (mixtures)</li>
              <li>• Lab conditions differ from field exposure</li>
            </ul>
          </div>
          <div className="bg-emerald-50 rounded-xl p-5 border border-emerald-200">
            <h4 className="font-bold text-emerald-800 mb-3">✅ Responsible Use</h4>
            <ul className="text-sm text-emerald-700 space-y-2">
              <li>• Use as screening tool, not final assessment</li>
              <li>• Apply precautionary principle when uncertain</li>
              <li>• Complement with experimental validation</li>
              <li>• Consider sub-lethal effects not captured here</li>
            </ul>
          </div>
        </div>
      </section>

      {/* References */}
      <section>
        <h3 className="text-xl font-bold text-journal-text mb-4">References & Further Reading</h3>
        <div className="bg-gray-50 rounded-xl p-5 border border-gray-200">
          <ul className="text-sm space-y-3">
            <li>
              <strong className="text-journal-text">ApisTox Dataset:</strong>{' '}
              <a href="https://www.nature.com/articles/s41597-024-04232-w" target="_blank" rel="noopener noreferrer" className="text-blue-600 hover:underline">
                Scientific Data (2024) — "ApisTox: A curated dataset for honey bee toxicity prediction"
              </a>
            </li>
            <li>
              <strong className="text-journal-text">SHAP Methodology:</strong>{' '}
              <a href="https://arxiv.org/abs/1705.07874" target="_blank" rel="noopener noreferrer" className="text-blue-600 hover:underline">
                Lundberg & Lee (2017) — "A Unified Approach to Interpreting Model Predictions"
              </a>
            </li>
            <li>
              <strong className="text-journal-text">XGBoost:</strong>{' '}
              <a href="https://arxiv.org/abs/1603.02754" target="_blank" rel="noopener noreferrer" className="text-blue-600 hover:underline">
                Chen & Guestrin (2016) — "XGBoost: A Scalable Tree Boosting System"
              </a>
            </li>
            <li>
              <strong className="text-journal-text">EPA Bee Toxicity Guidelines:</strong>{' '}
              <a href="https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/technical-overview-ecological-risk-assessment-0" target="_blank" rel="noopener noreferrer" className="text-blue-600 hover:underline">
                EPA Ecological Risk Assessment
              </a>
            </li>
          </ul>
        </div>
      </section>
    </div>
  );
};

// --- Main App Shell ---

const App = () => {
  const [activeTab, setActiveTab] = useState<TabView>('engine');

  return (
    <div className="min-h-screen bg-journal-bg text-journal-text font-sans selection:bg-journal-accent selection:text-white">
      {/* Navbar */}
      <nav className="border-b-2 border-journal-border bg-white/90 backdrop-blur-md sticky top-0 z-50 shadow-paper">
        <div className="max-w-7xl mx-auto px-6 h-20 flex items-center justify-between">
          <div className="flex items-center gap-4">
             <div className="w-12 h-12 bg-journal-accent/15 rounded-xl flex items-center justify-center text-journal-accent shadow-sm">
                <Flower2 size={28} />
             </div>
             <div>
                <h1 className="text-2xl font-display font-bold leading-none text-journal-text">ApisTox</h1>
                <p className="text-[10px] uppercase font-bold tracking-widest text-journal-dim">Pollinator Safety</p>
             </div>
          </div>

          <div className="flex gap-8">
            <button
              onClick={() => setActiveTab('engine')}
              className={`text-sm font-bold transition-all duration-quick border-b-3 py-6 ${activeTab === 'engine' ? 'border-journal-accent text-journal-accent' : 'border-transparent text-journal-dim hover:text-journal-text'}`}
            >
              Prediction Engine
            </button>
            <button
              onClick={() => setActiveTab('scenarios')}
              className={`text-sm font-bold transition-all duration-quick border-b-3 py-6 ${activeTab === 'scenarios' ? 'border-journal-accent text-journal-accent' : 'border-transparent text-journal-dim hover:text-journal-text'}`}
            >
              Scenario Analysis
            </button>
            <button
              onClick={() => setActiveTab('science')}
              className={`text-sm font-bold transition-all duration-quick border-b-3 py-6 ${activeTab === 'science' ? 'border-journal-accent text-journal-accent' : 'border-transparent text-journal-dim hover:text-journal-text'}`}
            >
              Scientific Basis
            </button>
          </div>
        </div>
      </nav>

      {/* Hero Section (only on home) */}
      {activeTab === 'engine' && (
        <header className="py-20 text-center border-b-2 border-journal-border/50 bg-journal-bg">
           <div className="max-w-4xl mx-auto px-6">
              <h1 className="text-5xl md:text-6xl font-display font-bold text-journal-text mb-6 leading-tight">
                 Pesticide Toxicity <br/>
                 <span className="italic text-journal-accent">Prediction for Apis mellifera</span>
              </h1>
              <p className="text-lg text-journal-dim font-serif max-w-2xl mx-auto leading-relaxed">
                 Machine learning classification of chemical compounds for honey bee toxicity assessment.
                 Built on the ApisTox dataset with XGBoost achieving 83.6% accuracy.
              </p>
           </div>
        </header>
      )}

      {/* Main Content */}
      <main className="max-w-7xl mx-auto px-6 py-12">
        {activeTab === 'engine' && <PredictionEngine />}
        {activeTab === 'scenarios' && <ScenariosView />}
        {activeTab === 'science' && <ScienceView />}
      </main>

      {/* Footer */}
      <footer className="border-t-2 border-journal-border mt-24 py-12 bg-white/70">
         <div className="max-w-7xl mx-auto px-6 text-center">
            <p className="font-serif text-2xl text-journal-text mb-4">ApisTox</p>
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
  );
};

// --- Custom Icons ---

const BugIcon = ({className}: {className?: string}) => (
  <svg className={className} width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="m8 2 1.88 1.88"/><path d="M14.12 3.88 16 2"/><path d="M9 7.13v-1a3.003 3.003 0 1 1 6 0v1"/><path d="M12 20c-3.3 0-6-2.7-6-6v-3a4 4 0 0 1 4-4h4a4 4 0 0 1 4 4v3c0 3.3-2.7 6-6 6"/><path d="M12 20v-9"/><path d="M6.53 9C4.6 8.8 3 7.1 3 5"/><path d="M6 13a6 6 0 0 0-6-6"/><path d="M18 13a6 6 0 0 1 6-6"/><path d="M17.47 9c1.93-.2 3.53-1.9 3.53-4"/></svg>
);

const BeeIcon = ({className}: {className?: string}) => (
  <svg className={className} width="24" height="24" viewBox="0 0 24 24" fill="currentColor" stroke="none"><path d="M19.4 8.7a5.05 5.05 0 0 0-1.2-1.2l-.7.7a3.96 3.96 0 0 1 1.2 1.2c.4.7.5 1.5.3 2.3l.9.3a5.03 5.03 0 0 0-.5-3.3z"/><path d="M4.6 8.7c.6-1.3 1.7-2.2 3-2.6l.3-.9a5.05 5.05 0 0 0-3.3.5 5.03 5.03 0 0 0-2 2.8l.9.3c.1-.8.5-1.5 1.1-2.1z"/><path d="M12 18c3.31 0 6-2.69 6-6s-2.69-6-6-6-6 2.69-6 6 2.69 6 6 6z" opacity=".3"/><path d="M12 4c-4.42 0-8 3.58-8 8s3.58 8 8 8 8-3.58 8-8-3.58-8-8-8zm0 14c-3.31 0-6-2.69-6-6s2.69-6 6-6 6 2.69 6 6-2.69 6-6 6z"/><path d="M12 6c-3.31 0-6 2.69-6 6s2.69 6 6 6 6-2.69 6-6-2.69-6-6-6zm0 10c-2.21 0-4-1.79-4-4s1.79-4 4-4 4 1.79 4 4-1.79 4-4 4z"/></svg>
);

export default App;
