import React, { useState } from 'react'

interface NavItem {
  id: string
  label: string
  icon: string
  path: string
}

interface SidebarProps {
  currentPath: string
  onNavigate: (path: string) => void
}

const navItems: NavItem[] = [
  { id: 'dashboard', label: 'Dashboard', icon: '📊', path: '/' },
  { id: 'explorer', label: 'Data Explorer', icon: '🔍', path: '/explorer' },
  { id: 'predict', label: 'Predict Toxicity', icon: '🎯', path: '/predict' },
  { id: 'model', label: 'Model Info', icon: '🤖', path: '/model' },
  { id: 'docs', label: 'Documentation', icon: '📚', path: '/docs' },
]

export const Sidebar: React.FC<SidebarProps> = ({ currentPath, onNavigate }) => {
  const [isOpen, setIsOpen] = useState(true)

  return (
    <>
      {/* Mobile toggle button */}
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="lg:hidden fixed top-4 left-4 z-50 p-2 rounded-lg bg-white shadow-md border border-border-default"
      >
        <span className="text-xl">{isOpen ? '✕' : '☰'}</span>
      </button>

      {/* Sidebar */}
      <aside
        className={`fixed left-0 top-0 h-full bg-white border-r border-border-default z-40 transition-transform duration-300 ${
          isOpen ? 'translate-x-0' : '-translate-x-full lg:translate-x-0'
        }`}
        style={{ width: '280px' }}
      >
        {/* Logo section */}
        <div className="p-6 border-b border-border-default">
          <div className="flex items-center gap-3">
            <span className="text-4xl">🐝</span>
            <div>
              <h2 className="text-xl font-bold text-text-primary">ApisTox</h2>
              <p className="text-xs text-text-tertiary">Bee Toxicity Predictor</p>
            </div>
          </div>
        </div>

        {/* Navigation */}
        <nav className="p-4">
          <ul className="space-y-1">
            {navItems.map((item) => {
              const isActive = currentPath === item.path
              return (
                <li key={item.id}>
                  <button
                    onClick={() => {
                      onNavigate(item.path)
                      // Close sidebar on mobile after navigation
                      if (window.innerWidth < 1024) {
                        setIsOpen(false)
                      }
                    }}
                    className={`w-full flex items-center gap-3 px-4 py-3 rounded-lg text-left transition-colors ${
                      isActive
                        ? 'bg-honey-light text-bee-black font-medium'
                        : 'text-text-secondary hover:bg-bg-secondary hover:text-text-primary'
                    }`}
                  >
                    <span className="text-xl">{item.icon}</span>
                    <span>{item.label}</span>
                  </button>
                </li>
              )
            })}
          </ul>
        </nav>

        {/* Footer info */}
        <div className="absolute bottom-0 left-0 right-0 p-4 border-t border-border-default">
          <div className="text-xs text-text-tertiary">
            <p className="font-medium mb-1">Model Performance</p>
            <div className="flex justify-between">
              <span>Accuracy:</span>
              <span className="font-semibold text-accent-safe">83.6%</span>
            </div>
            <div className="flex justify-between">
              <span>Dataset:</span>
              <span className="font-semibold">1,035</span>
            </div>
          </div>
        </div>
      </aside>

      {/* Overlay for mobile */}
      {isOpen && (
        <div
          className="lg:hidden fixed inset-0 bg-black/30 z-30"
          onClick={() => setIsOpen(false)}
        />
      )}
    </>
  )
}
