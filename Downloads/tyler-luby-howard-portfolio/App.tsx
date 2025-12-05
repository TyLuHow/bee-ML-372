import React, { useState } from 'react';
import { Navbar } from './components/Navbar';
import { Hero } from './components/Hero';
import { ProjectCard } from './components/ProjectCard';
import { PROJECTS, EXPERIENCE, RESOURCES, BLOG_POSTS, PAPERS } from './constants';
import { ArrowRight, Download, Mail, Linkedin, Github, Phone, FileText } from 'lucide-react';

function App() {
  const webProjects = PROJECTS.filter(p => p.category === 'Web App' || p.category === 'ML/AI');
  const enterpriseProjects = PROJECTS.filter(p => p.category === 'Enterprise');
  const [resourceFilter, setResourceFilter] = useState<'All' | 'Tool' | 'Tutorial' | 'Demo'>('All');

  const filteredResources = resourceFilter === 'All' 
    ? RESOURCES 
    : RESOURCES.filter(r => r.type === resourceFilter);

  return (
    <div className="min-h-screen bg-white font-sans text-slate-800">
      <Navbar />
      
      <main>
        <Hero />

        {/* --- PROJECTS SECTION --- */}
        <section id="projects" className="py-24 bg-white relative">
          <div className="absolute top-0 left-0 w-full h-px bg-gradient-to-r from-transparent via-slate-200 to-transparent"></div>
          
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="mb-16">
              <span className="font-mono text-fern text-sm tracking-widest uppercase mb-2 block">01 / Portfolio</span>
              <h2 className="text-3xl md:text-4xl font-bold text-slate-900">Selected Works</h2>
            </div>

            {/* Subsection 1: Web Apps */}
            <div className="mb-24">
              <div className="flex items-center mb-8 border-b border-slate-100 pb-4">
                <span className="font-mono text-slate-400 mr-4">A.</span>
                <h3 className="text-xl font-bold text-slate-800 tracking-tight">
                  Web Applications & AI Systems
                </h3>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-x-8 gap-y-12">
                {webProjects.map((project, idx) => (
                  <ProjectCard key={project.id} project={project} index={idx} />
                ))}
              </div>
            </div>

            {/* Subsection 2: Enterprise */}
            <div>
               <div className="flex items-center mb-8 border-b border-slate-100 pb-4">
                <span className="font-mono text-slate-400 mr-4">B.</span>
                <h3 className="text-xl font-bold text-slate-800 tracking-tight">
                  Enterprise Engineering
                </h3>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-x-8 gap-y-12">
                {enterpriseProjects.map((project, idx) => (
                  <ProjectCard key={project.id} project={project} index={idx} />
                ))}
              </div>
            </div>
          </div>
        </section>

        {/* --- EXPERIENCE SECTION --- */}
        <section id="experience" className="py-24 bg-slate-50/50 border-y border-slate-200">
          <div className="max-w-5xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="mb-16 flex flex-col md:flex-row md:items-end justify-between gap-4">
              <div>
                <span className="font-mono text-fern text-sm tracking-widest uppercase mb-2 block">02 / Career</span>
                <h2 className="text-3xl md:text-4xl font-bold text-slate-900">Professional Experience</h2>
              </div>
              <p className="text-slate-500 text-sm max-w-md text-right hidden md:block">
                Delivering quantified business impact through engineering and data science.
              </p>
            </div>

            <div className="relative pl-4 md:pl-0">
              {/* Technical Vertical Line */}
              <div className="absolute left-0 md:left-1/2 top-0 bottom-0 w-px bg-slate-300 border-l border-dashed border-slate-300 hidden md:block"></div>
              <div className="absolute left-0 top-0 bottom-0 w-px bg-slate-200 border-l border-dashed border-slate-300 md:hidden"></div>

              <div className="space-y-16">
                {EXPERIENCE.map((job, idx) => (
                  <div key={job.id} className={`relative flex flex-col md:flex-row gap-8 items-start ${idx % 2 === 0 ? 'md:flex-row-reverse' : ''}`}>
                    
                    {/* Center Node */}
                    <div className="absolute left-[-5px] md:left-1/2 md:-translate-x-[5px] w-[11px] h-[11px] bg-white border-2 border-fern rotate-45 z-10 mt-1.5"></div>

                    {/* Date (Opposite side on desktop) */}
                    <div className={`hidden md:block w-1/2 pt-1 ${idx % 2 === 0 ? 'text-left pl-12' : 'text-right pr-12'}`}>
                       <span className="font-mono text-sm text-slate-400 bg-white px-2 py-1 border border-slate-200 rounded-sm">
                         {job.period}
                       </span>
                    </div>

                    {/* Card Content */}
                    <div className="flex-grow w-full md:w-1/2 pl-6 md:pl-0">
                       <div className={`bg-white p-6 rounded-sm border border-slate-200 shadow-sm hover:border-fern/30 transition-all ${idx % 2 === 0 ? 'md:mr-12' : 'md:ml-12'}`}>
                         <div className="flex flex-col mb-4">
                            <div className="md:hidden font-mono text-xs text-slate-400 mb-1">{job.period}</div>
                            <h3 className="text-lg font-bold text-slate-900">{job.title}</h3>
                            <div className="text-fern font-medium text-sm font-mono mt-1 uppercase tracking-wide">{job.company}</div>
                         </div>
                         
                         <ul className="space-y-3 mb-5">
                           {job.achievements.map((ach, i) => (
                             <li key={i} className="text-slate-600 text-sm leading-relaxed flex items-start">
                               <span className="mr-2 text-slate-300 mt-1.5 text-[8px] font-mono">0{i+1}</span>
                               {ach}
                             </li>
                           ))}
                         </ul>

                         <div className="pt-4 border-t border-slate-100 border-dashed flex flex-wrap gap-2">
                           {job.techUsed.map(t => (
                             <span key={t} className="text-[10px] font-mono text-slate-500 bg-slate-50 px-2 py-1 rounded-sm border border-slate-100">
                               {t}
                             </span>
                           ))}
                         </div>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </section>

        {/* --- RESOURCES & BLOG --- */}
        <section id="resources" className="py-24 bg-white">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="grid grid-cols-1 lg:grid-cols-12 gap-12">
              
              {/* Left Col: Resources (7 cols) */}
              <div className="lg:col-span-7">
                <div className="flex flex-col sm:flex-row sm:items-end justify-between mb-10 gap-4">
                  <div>
                    <span className="font-mono text-fern text-sm tracking-widest uppercase mb-2 block">03 / Knowledge Base</span>
                    <h2 className="text-3xl font-bold text-slate-900">AI Resources Hub</h2>
                  </div>
                  {/* Technical Tabs */}
                  <div className="flex p-1 bg-slate-100 rounded-md">
                    {(['All', 'Tool', 'Tutorial', 'Demo'] as const).map(filter => (
                      <button
                        key={filter}
                        onClick={() => setResourceFilter(filter)}
                        className={`text-xs font-medium px-4 py-1.5 rounded-sm transition-all ${
                          resourceFilter === filter 
                            ? 'bg-white text-slate-900 shadow-sm' 
                            : 'text-slate-500 hover:text-slate-800'
                        }`}
                      >
                        {filter}
                      </button>
                    ))}
                  </div>
                </div>
                
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  {filteredResources.map(resource => (
                    <a key={resource.id} href={resource.link} className="block group h-full">
                      <div className="bg-white border border-slate-200 rounded-md p-5 hover:border-fern hover:shadow-md transition-all h-full flex flex-col">
                         <div className="flex justify-between items-start mb-3">
                            <span className={`text-[10px] font-mono uppercase tracking-wider px-2 py-0.5 rounded-sm border ${
                              resource.type === 'Tool' ? 'bg-blue-50 text-blue-600 border-blue-100' :
                              resource.type === 'Tutorial' ? 'bg-purple-50 text-purple-600 border-purple-100' :
                              'bg-orange-50 text-orange-600 border-orange-100'
                            }`}>
                                {resource.type}
                            </span>
                            <ArrowRight size={14} className="text-slate-300 group-hover:text-fern -translate-x-1 group-hover:translate-x-0 transition-all" />
                         </div>
                         <h4 className="font-bold text-slate-800 mb-2 group-hover:text-fern transition-colors">{resource.title}</h4>
                         <p className="text-sm text-slate-500 leading-relaxed flex-grow">{resource.description}</p>
                      </div>
                    </a>
                  ))}
                </div>
              </div>

              {/* Right Col: Blog & Papers (5 cols) */}
              <div className="lg:col-span-5 space-y-12">
                {/* Blog */}
                <div id="blog">
                   <div className="mb-8 border-b border-slate-100 pb-4">
                      <span className="font-mono text-fern text-sm tracking-widest uppercase mb-1 block">04 / Writing</span>
                      <h2 className="text-2xl font-bold text-slate-900">Recent Insights</h2>
                   </div>
                   <div className="space-y-6">
                     {BLOG_POSTS.map(post => (
                       <article key={post.id} className="group cursor-pointer border-l-2 border-transparent hover:border-fern pl-4 transition-all">
                         <div className="flex items-center text-[10px] font-mono text-slate-400 mb-1 space-x-2 uppercase">
                           <span>{post.date}</span>
                           <span>//</span>
                           <span>{post.readTime}</span>
                         </div>
                         <h3 className="text-base font-bold text-slate-800 mb-2 group-hover:text-fern transition-colors">
                           {post.title}
                         </h3>
                         <div className="flex gap-2">
                           {post.tags.map(tag => (
                             <span key={tag} className="text-[10px] text-slate-500 font-mono">
                               #{tag}
                             </span>
                           ))}
                         </div>
                       </article>
                     ))}
                   </div>
                </div>

                {/* Papers */}
                <div id="papers">
                  <div className="mb-6 border-b border-slate-100 pb-4">
                      <span className="font-mono text-fern text-sm tracking-widest uppercase mb-1 block">05 / Research</span>
                      <h2 className="text-2xl font-bold text-slate-900">Publications</h2>
                   </div>
                   <div className="space-y-4">
                       {PAPERS.map(paper => (
                          <div key={paper.id} className="bg-slate-50 p-4 rounded-md border border-slate-200 hover:border-fern/50 transition-colors group">
                             <div className="flex justify-between items-start">
                                <div>
                                   <h4 className="font-bold text-sm text-slate-800 mb-1 group-hover:text-fern transition-colors">{paper.title}</h4>
                                   <p className="text-xs text-slate-500 font-mono mb-2 uppercase">{paper.publication}</p>
                                </div>
                                <FileText size={14} className="text-slate-400 group-hover:text-fern" />
                             </div>
                          </div>
                       ))}
                    </div>
                </div>
              </div>
            </div>
          </div>
        </section>

        {/* --- RESUME SECTION --- */}
        <section id="resume" className="py-20 bg-slate-900 text-white relative overflow-hidden">
           {/* Background Pattern */}
           <div className="absolute inset-0 opacity-10" style={{ backgroundImage: 'radial-gradient(#ffffff 1px, transparent 1px)', backgroundSize: '24px 24px' }}></div>
           
           <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 relative z-10 text-center">
              <span className="font-mono text-fern text-sm tracking-widest uppercase mb-4 block">06 / Documentation</span>
              <h2 className="text-3xl md:text-4xl font-bold text-white mb-8">Curriculum Vitae</h2>
              
              <div className="bg-slate-800/50 backdrop-blur-sm border border-slate-700 rounded-lg p-8 md:p-12 max-w-2xl mx-auto hover:border-fern/50 transition-all duration-500 group">
                  <div className="flex flex-col items-center">
                     <div className="w-16 h-16 bg-slate-800 rounded-full flex items-center justify-center border border-slate-600 mb-6 group-hover:scale-110 transition-transform shadow-lg shadow-black/20">
                        <FileText size={28} className="text-fern" />
                     </div>
                     <h3 className="text-xl font-bold text-white mb-2">Tyler Luby Howard</h3>
                     <p className="text-slate-400 font-mono text-sm mb-8">Full-Stack Engineer & Industrial Operations</p>
                     
                     <div className="grid grid-cols-2 gap-x-12 gap-y-4 text-left text-sm text-slate-300 mb-8 border-t border-b border-slate-700/50 py-6 w-full">
                        <div>
                          <span className="block text-xs font-mono text-slate-500 uppercase mb-1">Experience</span>
                          3+ Years
                        </div>
                        <div>
                          <span className="block text-xs font-mono text-slate-500 uppercase mb-1">Education</span>
                          MS Industrial Eng.
                        </div>
                        <div>
                           <span className="block text-xs font-mono text-slate-500 uppercase mb-1">Focus</span>
                           AI Systems / Ops
                        </div>
                        <div>
                           <span className="block text-xs font-mono text-slate-500 uppercase mb-1">Location</span>
                           California, USA
                        </div>
                     </div>

                     <button className="group/btn relative inline-flex items-center justify-center px-8 py-3 font-bold text-white transition-all duration-200 bg-fern font-sans rounded-sm hover:bg-fern-light focus:outline-none ring-offset-2 focus:ring-2 ring-fern">
                        <Download size={18} className="mr-2 group-hover/btn:animate-bounce" />
                        Download Full Resume
                     </button>
                  </div>
              </div>
           </div>
        </section>

        {/* --- CONTACT SECTION --- */}
        <section id="contact" className="py-24 bg-white border-t border-slate-200">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="grid md:grid-cols-2 gap-16 items-center">
               <div>
                  <span className="font-mono text-fern text-sm tracking-widest uppercase mb-4 block">07 / Contact</span>
                  <h2 className="text-4xl md:text-5xl font-bold text-slate-900 mb-6 tracking-tight">Let's build something efficient.</h2>
                  <p className="text-lg text-slate-600 leading-relaxed mb-8 max-w-lg">
                    I am currently available for new opportunities in AI systems engineering and full-stack development.
                  </p>
                  
                  <a href="mailto:tyler@example.com" className="inline-flex items-center text-lg font-bold text-slate-900 hover:text-fern border-b-2 border-slate-200 hover:border-fern transition-all pb-1 mb-12">
                    tyler@example.com <ArrowRight size={20} className="ml-2" />
                  </a>
               </div>

               <div className="bg-slate-50 rounded-sm p-8 border border-slate-200">
                  <h3 className="font-mono text-xs font-bold text-slate-400 uppercase tracking-widest mb-8">Connect</h3>
                  <div className="space-y-6">
                     <a href="#" className="flex items-center justify-between group p-4 bg-white border border-slate-200 rounded-sm hover:border-fern hover:shadow-sm transition-all">
                        <div className="flex items-center">
                           <Linkedin size={20} className="text-slate-400 group-hover:text-[#0077b5] transition-colors mr-4" />
                           <span className="font-bold text-slate-700 group-hover:text-slate-900">LinkedIn</span>
                        </div>
                        <ArrowRight size={16} className="text-slate-300 group-hover:text-fern opacity-0 group-hover:opacity-100 transition-all" />
                     </a>
                     
                     <a href="#" className="flex items-center justify-between group p-4 bg-white border border-slate-200 rounded-sm hover:border-fern hover:shadow-sm transition-all">
                        <div className="flex items-center">
                           <Github size={20} className="text-slate-400 group-hover:text-slate-900 transition-colors mr-4" />
                           <span className="font-bold text-slate-700 group-hover:text-slate-900">GitHub</span>
                        </div>
                        <ArrowRight size={16} className="text-slate-300 group-hover:text-fern opacity-0 group-hover:opacity-100 transition-all" />
                     </a>

                     <a href="#" className="flex items-center justify-between group p-4 bg-white border border-slate-200 rounded-sm hover:border-fern hover:shadow-sm transition-all">
                        <div className="flex items-center">
                           <Phone size={20} className="text-slate-400 group-hover:text-fern transition-colors mr-4" />
                           <span className="font-bold text-slate-700 group-hover:text-slate-900">Schedule Call</span>
                        </div>
                        <ArrowRight size={16} className="text-slate-300 group-hover:text-fern opacity-0 group-hover:opacity-100 transition-all" />
                     </a>
                  </div>
               </div>
            </div>

            <div className="mt-24 pt-8 border-t border-slate-100 flex flex-col md:flex-row justify-between items-center text-xs text-slate-400 font-mono">
              <p>&copy; {new Date().getFullYear()} Tyler Luby Howard</p>
              <div className="flex space-x-6 mt-4 md:mt-0">
                 <span>San Luis Obispo, CA</span>
                 <span>Built with React & Tailwind</span>
              </div>
            </div>
          </div>
        </section>

      </main>
    </div>
  );
}

export default App;