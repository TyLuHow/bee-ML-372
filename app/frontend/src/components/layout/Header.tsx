import React from 'react'

interface HeaderProps {
  currentPath: string
}

const pathTitles: Record<string, { title: string; subtitle: string }> = {
  '/': {
    title: 'Dashboard',
    subtitle: 'Overview of pesticide toxicity predictions',
  },
  '/explorer': {
    title: 'Data Explorer',
    subtitle: 'Explore the ApisTox dataset',
  },
  '/predict': {
    title: 'Predict Toxicity',
    subtitle: 'ML-powered pesticide safety assessment',
  },
  '/model': {
    title: 'Model Information',
    subtitle: 'XGBoost classifier performance metrics',
  },
  '/docs': {
    title: 'Documentation',
    subtitle: 'User guide and API reference',
  },
}

export const Header: React.FC<HeaderProps> = ({ currentPath }) => {
  const { title, subtitle } = pathTitles[currentPath] || pathTitles['/']

  return (
    <header className="bg-white border-b border-border-default px-8 py-6">
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold text-text-primary">{title}</h1>
          <p className="text-sm text-text-secondary mt-1">{subtitle}</p>
        </div>

        <div className="flex items-center gap-4">
          {/* Future: Add user menu, settings, etc */}
          <div className="flex items-center gap-2 px-4 py-2 bg-bg-secondary rounded-lg">
            <span className="text-sm text-text-secondary">IME 372 Project</span>
          </div>
        </div>
      </div>
    </header>
  )
}
