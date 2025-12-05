import React, { useState } from 'react';
import { Project } from '../types';
import { ExternalLink, X, ArrowUpRight } from 'lucide-react';

interface ProjectCardProps {
  project: Project;
  index: number;
}

export const ProjectCard: React.FC<ProjectCardProps> = ({ project, index }) => {
  const [isModalOpen, setIsModalOpen] = useState(false);
  const projectCode = `${project.category === 'Web App' ? 'WEB' : project.category === 'Enterprise' ? 'ENT' : 'ML'} // ${String(index + 1).padStart(2, '0')}`;

  return (
    <>
      <div 
        className="group bg-white rounded-md border border-slate-200 overflow-hidden hover:border-fern/50 hover:shadow-lg transition-all duration-300 cursor-pointer flex flex-col h-full relative"
        onClick={() => setIsModalOpen(true)}
      >
        {/* Technical Header Bar */}
        <div className="bg-slate-50 border-b border-slate-100 px-4 py-2 flex justify-between items-center">
          <span className="font-mono text-[10px] text-slate-500 tracking-wider uppercase">{projectCode}</span>
          <div className="flex space-x-1">
            <div className="w-2 h-2 rounded-full bg-slate-200 group-hover:bg-fern/40 transition-colors"></div>
            <div className="w-2 h-2 rounded-full bg-slate-200 group-hover:bg-fern/40 transition-colors delay-75"></div>
          </div>
        </div>

        {/* Card Image */}
        <div className="relative h-48 overflow-hidden bg-slate-100 border-b border-slate-100 group">
          <img 
            src={project.image} 
            alt={project.title} 
            className="w-full h-full object-cover grayscale-[20%] group-hover:grayscale-0 transition-all duration-500 group-hover:scale-105"
          />
          <div className="absolute top-2 right-2 opacity-0 group-hover:opacity-100 transition-opacity duration-300">
             <div className="bg-white/90 p-1.5 rounded-sm shadow-sm">
               <ArrowUpRight size={16} className="text-fern" />
             </div>
          </div>
        </div>

        {/* Card Body */}
        <div className="p-5 flex flex-col flex-grow relative">
          <h3 className="text-lg font-bold text-slate-900 mb-2 leading-tight group-hover:text-fern transition-colors">{project.title}</h3>
          <p className="text-slate-600 text-sm mb-4 line-clamp-3 leading-relaxed flex-grow">
            {project.shortDescription}
          </p>
          
          <div className="mt-auto pt-4 border-t border-slate-100 border-dashed">
            {/* Metric */}
            <div className="flex items-center text-xs font-medium text-slate-700 mb-2">
              <span className="w-1.5 h-1.5 bg-fern rounded-sm mr-2"></span>
              {project.metrics[0]}
            </div>
            {/* Stack */}
            <div className="text-[10px] font-mono text-slate-400 truncate uppercase tracking-tight">
              {project.techStack.join(' / ')}
            </div>
          </div>
        </div>
      </div>

      {/* Modal Overlay */}
      {isModalOpen && (
        <div className="fixed inset-0 z-[100] flex items-center justify-center px-4 py-6 sm:px-6">
          <div 
            className="absolute inset-0 bg-slate-900/80 backdrop-blur-sm transition-opacity"
            onClick={() => setIsModalOpen(false)}
          />
          
          <div className="relative bg-white rounded-lg shadow-2xl max-w-5xl w-full max-h-[90vh] overflow-y-auto flex flex-col md:flex-row overflow-hidden animate-[fadeIn_0.2s_ease-out] border border-slate-200">
            <button 
              onClick={(e) => { e.stopPropagation(); setIsModalOpen(false); }}
              className="absolute top-4 right-4 z-20 p-2 bg-white/90 rounded-full text-slate-600 hover:text-fern hover:rotate-90 transition-all shadow-sm"
            >
              <X size={20} />
            </button>

            {/* Modal Image/Demo */}
            <div className="w-full md:w-3/5 bg-slate-100 relative min-h-[300px] md:min-h-full border-b md:border-b-0 md:border-r border-slate-200">
                 <img src={project.image} alt={project.title} className="w-full h-full object-cover" />
                 {project.link && (
                   <a 
                    href={project.link}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="absolute bottom-6 right-6 bg-white/95 text-slate-900 px-5 py-2.5 rounded-sm text-xs font-bold uppercase tracking-wide shadow-lg hover:bg-fern hover:text-white transition-all flex items-center gap-2 border border-slate-200"
                   >
                     Live Project <ExternalLink size={12} />
                   </a>
                 )}
            </div>
            
            {/* Modal Content */}
            <div className="w-full md:w-2/5 p-8 md:p-10 flex flex-col h-full bg-white overflow-y-auto">
                <div className="mb-6">
                  <span className="font-mono text-xs text-fern uppercase tracking-wider mb-2 block border-l-2 border-fern pl-3">
                    {project.category}
                  </span>
                  <h2 className="text-3xl font-bold text-slate-900 tracking-tight">{project.title}</h2>
                </div>

                <div className="space-y-8 flex-grow">
                  <div>
                    <h4 className="font-mono text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Context</h4>
                    <p className="text-slate-600 text-sm leading-relaxed">
                      {project.fullDescription}
                    </p>
                  </div>
                  
                  <div>
                    <h4 className="font-mono text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Impact</h4>
                    <ul className="space-y-3">
                      {project.metrics.map((m, i) => (
                        <li key={i} className="flex items-start text-sm text-slate-800 font-medium">
                          <span className="mr-3 text-fern mt-0.5 text-xs font-mono">{`0${i+1}`}</span>
                          {m}
                        </li>
                      ))}
                    </ul>
                  </div>

                  <div>
                    <h4 className="font-mono text-xs font-bold text-slate-400 uppercase tracking-widest mb-3">Technology</h4>
                    <div className="flex flex-wrap gap-2">
                      {project.techStack.map((tech) => (
                        <span key={tech} className="px-2 py-1 bg-slate-50 border border-slate-200 rounded-sm text-[10px] font-mono text-slate-600 uppercase">
                          {tech}
                        </span>
                      ))}
                    </div>
                  </div>
                </div>
            </div>
          </div>
        </div>
      )}
    </>
  );
};