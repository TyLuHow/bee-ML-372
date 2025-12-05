import React, { useState, useEffect } from 'react';
import { ArrowRight, ArrowLeft } from 'lucide-react';
import { PROJECTS } from '../constants';

// Filter only the featured web apps for the hero
const FEATURED_PROJECTS = PROJECTS.filter(p => ['stormwater', 'onroute', 'hypertrack', 'beeml'].includes(p.id));

export const Hero: React.FC = () => {
  const [currentSlide, setCurrentSlide] = useState(0);

  useEffect(() => {
    const timer = setInterval(() => {
      setCurrentSlide((prev) => (prev + 1) % FEATURED_PROJECTS.length);
    }, 6000);
    return () => clearInterval(timer);
  }, []);

  const nextSlide = () => setCurrentSlide((prev) => (prev + 1) % FEATURED_PROJECTS.length);
  const prevSlide = () => setCurrentSlide((prev) => (prev - 1 + FEATURED_PROJECTS.length) % FEATURED_PROJECTS.length);

  return (
    <section id="hero" className="relative pt-24 pb-16 md:pt-32 md:pb-24 overflow-hidden bg-white">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        
        {/* Intro Text */}
        <div className="text-center mb-16">
          <h1 className="text-5xl md:text-7xl font-bold text-slate-900 mb-6 tracking-tight">
            Tyler Luby Howard
          </h1>
          <p className="text-lg md:text-xl text-slate-500 font-light max-w-2xl mx-auto leading-relaxed">
            <span className="text-slate-900 font-medium">Industrial Engineer</span> meeting <span className="text-slate-900 font-medium">Full-Stack Architecture</span>.
            <br className="hidden md:block"/> Building scalable AI systems and data-driven enterprise solutions.
          </p>
        </div>

        {/* Carousel Container */}
        <div className="relative w-full max-w-6xl mx-auto aspect-[16/9] md:aspect-[21/9] bg-white rounded-md overflow-hidden shadow-xl border border-slate-200 group">
          {FEATURED_PROJECTS.map((project, index) => (
            <div
              key={project.id}
              className={`absolute inset-0 w-full h-full transition-opacity duration-700 ease-in-out ${
                index === currentSlide ? 'opacity-100 z-10' : 'opacity-0 z-0'
              }`}
            >
              {/* Content Layout */}
              <div className="flex flex-col md:flex-row h-full">
                
                {/* Image Section */}
                <div className="w-full md:w-7/12 h-2/3 md:h-full relative overflow-hidden bg-slate-100 border-r border-slate-100">
                  <img 
                    src={project.image} 
                    alt={project.title} 
                    className="object-cover w-full h-full transform transition-transform duration-[10000ms] scale-100 group-hover:scale-105"
                  />
                </div>

                {/* Text/Info Section */}
                <div className="w-full md:w-5/12 h-1/3 md:h-full bg-white p-8 md:p-12 flex flex-col justify-center relative z-20">
                  <div className="mb-auto hidden md:block">
                     <span className="font-mono text-xs text-slate-400 uppercase tracking-widest">Featured Project 0{index + 1}</span>
                  </div>
                  
                  <div>
                    <h3 className="text-2xl md:text-3xl font-bold text-slate-900 mb-3">
                      {project.title}
                    </h3>
                    <p className="text-sm md:text-base text-slate-600 mb-6 leading-relaxed">
                      {project.shortDescription}
                    </p>
                    <div className="flex flex-wrap gap-2 mb-8">
                       {project.techStack.slice(0,3).map(tech => (
                         <span key={tech} className="font-mono text-[10px] bg-slate-50 border border-slate-100 px-2 py-1 rounded text-slate-500 uppercase">{tech}</span>
                       ))}
                    </div>
                  </div>
                  
                  <div className="mt-auto">
                    <a 
                      href="#projects"
                      className="inline-flex items-center text-sm font-bold text-white bg-slate-900 px-6 py-3 rounded-sm hover:bg-fern transition-colors group/btn"
                    >
                      View Details
                      <ArrowRight size={16} className="ml-2 transform group-hover/btn:translate-x-1 transition-transform" />
                    </a>
                  </div>
                </div>
              </div>
            </div>
          ))}

          {/* Controls - Minimal Tech Style */}
          <div className="absolute bottom-6 left-6 md:left-auto md:right-12 z-30 flex space-x-3">
            {FEATURED_PROJECTS.map((_, idx) => (
              <button
                key={idx}
                onClick={() => setCurrentSlide(idx)}
                className={`h-1 transition-all duration-300 ${
                  idx === currentSlide ? 'bg-fern w-8' : 'bg-slate-200 w-4 hover:bg-slate-300'
                }`}
                aria-label={`Go to slide ${idx + 1}`}
              />
            ))}
          </div>

          <button 
            onClick={prevSlide}
            className="absolute left-4 top-1/2 -translate-y-1/2 z-30 p-3 rounded-full bg-white/90 hover:bg-white text-slate-700 shadow-sm border border-slate-100 md:hidden"
          >
            <ArrowLeft size={20} />
          </button>
          <button 
            onClick={nextSlide}
            className="absolute right-4 top-1/2 -translate-y-1/2 z-30 p-3 rounded-full bg-white/90 hover:bg-white text-slate-700 shadow-sm border border-slate-100 md:hidden"
          >
            <ArrowRight size={20} />
          </button>

        </div>
      </div>
    </section>
  );
};