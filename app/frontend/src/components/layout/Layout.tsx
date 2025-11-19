import React from 'react'
import { Sidebar } from './Sidebar'
import { Header } from './Header'

interface LayoutProps {
  children: React.ReactNode
  currentPath: string
  onNavigate: (path: string) => void
}

export const Layout: React.FC<LayoutProps> = ({ children, currentPath, onNavigate }) => {
  return (
    <div className="min-h-screen bg-bg-secondary">
      {/* Sidebar */}
      <Sidebar currentPath={currentPath} onNavigate={onNavigate} />

      {/* Main content area */}
      <div className="lg:ml-[280px]">
        {/* Header */}
        <Header currentPath={currentPath} />

        {/* Page content */}
        <main className="p-8">
          {children}
        </main>

        {/* Footer */}
        <footer className="border-t border-border-default px-8 py-6 bg-white">
          <div className="text-center text-text-tertiary text-sm">
            <p>Built for pollinator conservation</p>
            <p className="mt-1">
              Data: ApisTox Dataset (1,035 compounds) | Model: XGBoost | Interpretability: SHAP
            </p>
          </div>
        </footer>
      </div>
    </div>
  )
}
