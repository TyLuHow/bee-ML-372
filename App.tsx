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
  Play
} from 'lucide-react';
import { ChemicalData, PredictionResult, TabView, Scenario } from './types';
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
  const [formData, setFormData] = useState<ChemicalData>({
    name: '',
    mw: 350.5,
    logP: 3.2,
    exposure: 'Contact (Direct Spray)',
    category: 'Insecticide'
  });

  const handleAnalyze = async () => {
    setLoading(true);
    const res = await analyzeChemicalToxicity(formData);
    setResult(res);
    setLoading(false);
  };

  return (
    <div className="animate-fade-in">
      <SectionHeader 
        title="Analyze a Compound" 
        subtitle="Input chemical properties below to generate a real-time toxicity risk assessment." 
      />

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-8 items-start">
        {/* Input Form */}
        <div className="lg:col-span-5">
          <Card className="space-y-6">
            <div className="border-b border-journal-border pb-4 mb-4">
              <h3 className="font-bold text-xs uppercase tracking-widest text-journal-dim">1. Chemical Identity</h3>
            </div>
            
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-bold text-journal-text mb-2">Compound Name</label>
                <div className="relative">
                  <input
                    type="text"
                    placeholder="e.g. Imidacloprid"
                    value={formData.name}
                    onChange={(e) => setFormData({...formData, name: e.target.value})}
                    className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg focus:border-journal-accent focus:bg-white focus:shadow-sm outline-none transition-all duration-quick font-medium"
                  />
                  <Search size={18} className="absolute right-4 top-3.5 text-journal-dim/40" />
                </div>
              </div>

              <div className="grid grid-cols-2 gap-4">
                <div>
                   <label className="block text-sm font-bold text-journal-text mb-2">Data Source</label>
                   <select className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg outline-none text-sm font-medium focus:border-journal-accent focus:bg-white focus:shadow-sm transition-all duration-quick cursor-pointer">
                     <option>PPDB (Pesticide Properties)</option>
                     <option>ECOTOX (EPA)</option>
                   </select>
                </div>
                 <div>
                   <label className="block text-sm font-bold text-journal-text mb-2">Exposure Route</label>
                   <select
                      value={formData.exposure}
                      onChange={(e) => setFormData({...formData, exposure: e.target.value})}
                      className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg outline-none text-sm font-medium focus:border-journal-accent focus:bg-white focus:shadow-sm transition-all duration-quick cursor-pointer"
                   >
                     <option>Contact (Direct Spray)</option>
                     <option>Oral (Ingestion)</option>
                     <option>Systemic</option>
                   </select>
                </div>
              </div>
            </div>

            <div className="border-b border-journal-border pb-4 pt-2 mb-4">
              <h3 className="font-bold text-xs uppercase tracking-widest text-journal-dim">2. Pesticide Class</h3>
            </div>

            <div className="grid grid-cols-2 gap-3">
              {['Insecticide', 'Herbicide', 'Fungicide', 'Other'].map(cat => (
                <button
                  key={cat}
                  onClick={() => setFormData({...formData, category: cat})}
                  className={`p-5 border-2 rounded-xl flex flex-col items-center gap-2 transition-all duration-quick ${
                    formData.category === cat
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

            <div className="border-b border-journal-border pb-4 pt-2 mb-4">
              <h3 className="font-bold text-xs uppercase tracking-widest text-journal-dim">3. Molecular Fingerprint</h3>
            </div>

            <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-bold text-journal-text mb-2">Molecular Weight (g/mol)</label>
                  <input
                    type="number"
                    value={formData.mw}
                    onChange={(e) => setFormData({...formData, mw: parseFloat(e.target.value)})}
                    className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg focus:border-journal-accent focus:bg-white focus:shadow-sm outline-none transition-all duration-quick font-mono text-lg"
                  />
                  <span className="text-[10px] text-journal-dim mt-1.5 block leading-tight font-serif italic">Heavier molecules may have lower bioavailability.</span>
                </div>
                <div>
                  <label className="block text-sm font-bold text-journal-text mb-2">LogP (Lipophilicity)</label>
                  <input
                    type="number"
                    value={formData.logP}
                    onChange={(e) => setFormData({...formData, logP: parseFloat(e.target.value)})}
                    className="w-full px-4 py-3 bg-[#FFFAF0] border-b-2 border-journal-border rounded-lg focus:border-journal-accent focus:bg-white focus:shadow-sm outline-none transition-all duration-quick font-mono text-lg"
                  />
                  <span className="text-[10px] text-journal-dim mt-1.5 block leading-tight font-serif italic">High LogP often means higher toxicity accumulation.</span>
                </div>
            </div>

            <div className="pt-4">
              <Button onClick={handleAnalyze} disabled={loading} className="w-full text-lg">
                {loading ? <Activity className="animate-spin" /> : 'Run Prediction Model'}
              </Button>
            </div>
          </Card>
        </div>

        {/* Output */}
        <div className="lg:col-span-7 h-full">
           <Card className={`h-full min-h-[600px] flex flex-col justify-center items-center text-center transition-all duration-gentle ${!result ? 'bg-journal-surface' : (result.toxicity === 'Safe' ? 'bg-[#F0FDF4] border-[#86EFAC]' : 'bg-[#FFF1F2] border-[#FECACA]')}`}>
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
                      <div className={`flex flex-col items-center ${result.toxicity === 'Safe' ? 'text-journal-green' : 'text-journal-red'}`}>
                        {result.toxicity === 'Safe' ? <Leaf size={64} strokeWidth={1.5} /> : <AlertTriangle size={64} strokeWidth={1.5} />}
                        <h2 className="text-5xl font-display font-bold mt-4">{result.toxicity === 'Safe' ? 'Likely Safe' : 'Toxic'}</h2>
                        <p className="font-serif italic text-journal-dim mt-2">Molecular profile suggests {result.toxicity === 'Safe' ? 'low' : 'high'} risk to honey bees.</p>
                      </div>
                   </div>

                   <div className="w-full bg-white/70 rounded-xl p-6 mb-6 border border-black/5 shadow-sm">
                      <div className="flex justify-between items-end mb-3">
                         <span className="text-xs font-bold uppercase tracking-widest text-journal-dim">Model Confidence</span>
                         <span className="font-mono text-3xl font-bold">{result.confidence}%</span>
                      </div>
                      <div className="w-full h-3 bg-gray-200/80 rounded-full overflow-hidden shadow-inner">
                        <div
                          className={`h-full transition-all duration-gentle ${result.toxicity === 'Safe' ? 'bg-gradient-to-r from-journal-green to-emerald-500' : 'bg-gradient-to-r from-journal-red to-rose-600'}`}
                          style={{ width: `${result.confidence}%` }}
                        ></div>
                      </div>
                      <div className="flex justify-between mt-2 text-[10px] text-gray-400 font-mono uppercase tracking-wider">
                        <span>Uncertain</span>
                        <span>Certain</span>
                      </div>
                   </div>

                   <div className="space-y-6 w-full">
                     <div>
                       <div className="flex items-center gap-3 mb-3 text-journal-accent">
                          <Microscope size={22} />
                          <h4 className="font-bold font-display text-xl">Analysis</h4>
                       </div>
                       <p className="text-journal-text leading-relaxed pl-9">{result.explanation}</p>
                     </div>

                     <div>
                       <div className="flex items-center gap-3 mb-3 text-journal-accent">
                          <BookOpen size={22} />
                          <h4 className="font-bold font-display text-xl">Recommendation</h4>
                       </div>
                       <p className="text-journal-text leading-relaxed pl-9">{result.recommendation}</p>
                     </div>
                   </div>

                   <div className="mt-auto w-full pt-8 text-right">
                      <span className="font-mono text-xs text-journal-dim/50 uppercase">Model: XGBoost v1.4 | {new Date().toLocaleTimeString()}</span>
                   </div>
                </div>
              )}
           </Card>
        </div>
      </div>
    </div>
  );
};

const ScenariosView = () => {
  const [activeCandidate, setActiveCandidate] = useState<'A' | 'B' | null>(null);

  const scenarios: Scenario[] = [
    {
      id: 'almond',
      title: "The Almond Grower's Decision",
      role: "Maria Rodriguez",
      description: "Balancing pest control with pollinator safety during bloom season.",
      icon: "🌸"
    },
    {
      id: 'startup',
      title: "The R&D Startup's $3M Decision",
      role: "NexGen AgroSciences",
      description: "Prioritizing fungicide candidates before committing to expensive EPA testing.",
      icon: "⚗️"
    },
    {
      id: 'beekeeper',
      title: "The Beekeeper's Investigation",
      role: "James Patterson",
      description: "Analyzing hive collapse potential near treated sunflower fields.",
      icon: "🔍"
    }
  ];

  return (
    <div className="animate-fade-in space-y-12">
      <SectionHeader 
        title="Real-World Scenarios" 
        subtitle="Explore how ApisTox helps different stakeholders make data-driven decisions about pesticide safety." 
      />

      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        {scenarios.map((s) => (
          <div key={s.id} className={`p-8 border-2 rounded-xl bg-white shadow-card-subtle hover:shadow-lifted hover:border-journal-accent/40 transition-all duration-smooth cursor-pointer group ${s.id === 'startup' ? 'border-journal-accent ring-2 ring-journal-accent/20' : 'border-journal-border'}`}>
             <div className="text-5xl mb-5 group-hover:scale-110 transition-transform duration-smooth">{s.icon}</div>
             <h3 className="font-bold font-display text-xl mb-2 leading-tight">{s.title}</h3>
             <p className="text-xs font-bold text-journal-dim uppercase tracking-wider mb-4">{s.role}</p>
             <p className="text-sm text-journal-text leading-relaxed">{s.description}</p>
          </div>
        ))}
      </div>

      <div className="bg-[#FFFBF0] border border-journal-accent/20 rounded-2xl p-8 md:p-12 relative overflow-hidden">
        {/* Decorative Background Icon */}
        <div className="absolute top-0 right-0 opacity-5 pointer-events-none">
          <Beaker size={400} />
        </div>

        <div className="relative z-10">
          <div className="flex items-center gap-4 mb-3">
            <span className="text-5xl">⚗️</span>
            <h2 className="text-3xl font-display font-bold text-journal-text">The R&D Startup's $3M Decision</h2>
          </div>
          <p className="text-journal-maroon font-bold mb-8 text-lg">NexGen AgroSciences, Fungicide Development Team</p>
          
          <div className="max-w-3xl space-y-4 mb-8">
            <p className="text-lg leading-relaxed">
              The team has synthesized 40 novel fungicide candidates. Before committing <span className="font-bold">$62K per compound</span> for EPA bee toxicity testing (OECD 213/214/245), they need to prioritize which candidates to advance.
            </p>
            <div className="inline-flex items-center gap-2 bg-journal-accent/10 text-journal-accent px-4 py-2 rounded-full font-bold text-sm">
              <span>💰 $2.5M potential savings + 6-month timeline acceleration</span>
            </div>
          </div>

          <div className="flex flex-wrap gap-4 mb-8">
            <button
              onClick={() => setActiveCandidate('A')}
              className={`px-6 py-3 rounded-xl font-bold transition-all duration-quick ${activeCandidate === 'A' ? 'bg-journal-accent text-white shadow-lifted' : 'bg-white border-2 border-journal-border hover:bg-[#FFFAF0] hover:border-journal-accent/50'}`}
            >
              Candidate A (Triazole-based)
            </button>
            <button
              onClick={() => setActiveCandidate('B')}
              className={`px-6 py-3 rounded-xl font-bold transition-all duration-quick ${activeCandidate === 'B' ? 'bg-journal-accent text-white shadow-lifted' : 'bg-white border-2 border-journal-border hover:bg-[#FFFAF0] hover:border-journal-accent/50'}`}
            >
              Candidate B (Pyrimidine-based)
            </button>
            <button className="px-6 py-3 rounded-xl bg-journal-maroon text-white font-bold hover:bg-[#5C1D0A] transition-all duration-quick flex items-center gap-2 ml-auto shadow-card-subtle hover:shadow-lifted">
              Run All <Play size={16} fill="currentColor" />
            </button>
          </div>

          {activeCandidate && (
            <div className="bg-white rounded-xl p-8 border-2 border-journal-border shadow-lifted animate-fade-in flex gap-8">
              <div className="flex-1">
                <h4 className="font-bold text-2xl mb-3 font-display">
                  {activeCandidate === 'A' ? 'Candidate A (Triazole-based)' : 'Candidate B (Pyrimidine-based)'}
                </h4>
                <p className="text-journal-dim mb-5 leading-relaxed">
                  {activeCandidate === 'A'
                    ? 'High LogP (4.2), multiple aromatic rings - structural alerts present.'
                    : 'Low LogP (1.8), novel fused ring system.'}
                </p>
                <div className="flex items-center gap-3">
                   <span className="text-xs font-bold uppercase text-journal-dim tracking-wider">Expected:</span>
                   <span className={`px-3 py-1.5 rounded-lg text-xs font-bold ${activeCandidate === 'A' ? 'bg-red-100 text-red-700' : 'bg-green-100 text-green-700'}`}>
                     {activeCandidate === 'A' ? 'TOXIC' : 'LIKELY SAFE'}
                   </span>
                </div>
              </div>
              <div className="w-1/3 border-l-2 pl-8 border-gray-100 flex flex-col justify-center">
                 <div className="text-sm font-bold text-journal-dim mb-2">Confidence</div>
                 <div className="text-4xl font-mono font-bold mb-3">{activeCandidate === 'A' ? '87%' : '92%'}</div>
                 <div className="w-full bg-gray-100 h-3 rounded-full overflow-hidden shadow-inner">
                    <div className={`h-full transition-all duration-gentle ${activeCandidate === 'A' ? 'bg-gradient-to-r from-red-500 to-rose-600' : 'bg-gradient-to-r from-green-500 to-emerald-600'}`} style={{width: activeCandidate === 'A' ? '87%' : '92%'}}></div>
                 </div>
              </div>
            </div>
          )}

        </div>
      </div>
    </div>
  );
};

const ScienceView = () => {
  return (
    <div className="animate-fade-in space-y-16">
      <div className="text-center max-w-2xl mx-auto mb-12">
         <h2 className="text-4xl font-display font-bold text-journal-text mb-5">Understanding the Data</h2>
         <p className="text-lg text-journal-dim font-serif italic leading-relaxed">
           Our model learns from real-world toxicity data collected from three major scientific databases. Here's what powers our predictions.
         </p>
      </div>

      <div>
        <div className="flex items-center gap-4 mb-10">
           <span className="w-10 h-10 rounded-full bg-journal-accent/20 text-journal-accent font-bold flex items-center justify-center font-mono">1</span>
           <h3 className="text-2xl font-bold font-display">What We're Predicting</h3>
        </div>

        <div className="grid md:grid-cols-2 gap-8 bg-[#FAFAF8] p-10 rounded-2xl border-2 border-gray-100">
           <div className="bg-white p-10 rounded-xl text-center shadow-card-subtle border-2 border-emerald-100">
              <div className="inline-flex p-4 rounded-full bg-emerald-100 text-emerald-600 mb-5">
                 <CheckCircle2 size={36} />
              </div>
              <h4 className="text-2xl font-bold font-display text-emerald-800 mb-2">Non-Toxic</h4>
              <p className="font-mono text-emerald-600 mb-3 text-lg">LD50 &gt; 11 μg/bee</p>
              <p className="text-sm text-gray-500 font-serif italic">Safe for pollinators</p>
           </div>
           <div className="bg-white p-10 rounded-xl text-center shadow-card-subtle border-2 border-rose-100">
              <div className="inline-flex p-4 rounded-full bg-rose-100 text-rose-600 mb-5">
                 <AlertTriangle size={36} />
              </div>
              <h4 className="text-2xl font-bold font-display text-rose-800 mb-2">Toxic</h4>
              <p className="font-mono text-rose-600 mb-3 text-lg">LD50 ≤ 11 μg/bee</p>
              <p className="text-sm text-gray-500 font-serif italic">Harmful to honey bees</p>
           </div>
           <div className="md:col-span-2 text-center pt-4">
              <p className="text-sm font-bold text-journal-dim">LD50 = Lethal Dose that kills 50% of test population (EPA standard threshold)</p>
           </div>
        </div>
      </div>

      <div>
        <div className="flex items-center gap-4 mb-10">
           <span className="w-10 h-10 rounded-full bg-journal-accent/20 text-journal-accent font-bold flex items-center justify-center font-mono">2</span>
           <h3 className="text-2xl font-bold font-display">Data Sources</h3>
        </div>

        <div className="grid md:grid-cols-3 gap-6">
           <div className="p-8 bg-blue-50 rounded-xl border-2 border-blue-100 shadow-card-subtle">
              <h4 className="font-bold font-display text-blue-900 mb-2 text-xl">ECOTOX</h4>
              <p className="text-xs font-bold uppercase text-blue-500 mb-4 tracking-wider">EPA Ecotoxicology Database</p>
              <p className="text-sm text-blue-800 leading-relaxed">
                U.S. Environmental Protection Agency's curated database of chemical toxicity studies on aquatic and terrestrial life.
              </p>
           </div>
           <div className="p-8 bg-green-50 rounded-xl border-2 border-green-100 shadow-card-subtle">
              <h4 className="font-bold font-display text-green-900 mb-2 text-xl">PPDB</h4>
              <p className="text-xs font-bold uppercase text-green-600 mb-4 tracking-wider">Pesticide Properties Database</p>
              <p className="text-sm text-green-800 leading-relaxed">
                 University of Hertfordshire's comprehensive database covering physicochemical and ecotoxicological data.
              </p>
           </div>
           <div className="p-8 bg-purple-50 rounded-xl border-2 border-purple-100 shadow-card-subtle">
              <h4 className="font-bold font-display text-purple-900 mb-2 text-xl">BPDB</h4>
              <p className="text-xs font-bold uppercase text-purple-500 mb-4 tracking-wider">Bio-Pesticides Database</p>
              <p className="text-sm text-purple-800 leading-relaxed">
                 Specialized database for biological pesticides and their environmental impact profiles.
              </p>
           </div>
        </div>
      </div>
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
         <div className="max-w-7xl mx-auto px-6 flex flex-col md:flex-row justify-between items-center opacity-70">
            <div className="flex items-center gap-3 mb-4 md:mb-0">
               <Leaf size={18} className="text-journal-green" />
               <span className="font-bold text-sm font-display">ApisTox Project</span>
            </div>
            <p className="font-serif italic text-sm">Design inspired by nature. Powered by science.</p>
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