import { useState, useEffect, lazy, Suspense } from 'react'
import { Layout } from './components/layout/Layout'
import { AppProvider } from './store/AppContext'
import './App.css'

// Code splitting: Lazy load pages to reduce initial bundle size
const Dashboard = lazy(() => import('./pages/Dashboard').then(m => ({ default: m.Dashboard })))
const ExplorerPage = lazy(() => import('./pages/ExplorerPage').then(m => ({ default: m.ExplorerPage })))
const PredictPage = lazy(() => import('./pages/PredictPage').then(m => ({ default: m.PredictPage })))
const ModelPage = lazy(() => import('./pages/ModelPage').then(m => ({ default: m.ModelPage })))
const DocsPage = lazy(() => import('./pages/DocsPage').then(m => ({ default: m.DocsPage })))

function App() {
  const [currentPath, setCurrentPath] = useState('/')

  // Simple hash-based router
  useEffect(() => {
    const handleHashChange = () => {
      const hash = window.location.hash.slice(1) || '/'
      setCurrentPath(hash)
    }

    // Set initial path
    handleHashChange()

    // Listen for hash changes
    window.addEventListener('hashchange', handleHashChange)

    return () => {
      window.removeEventListener('hashchange', handleHashChange)
    }
  }, [])

  const navigate = (path: string) => {
    window.location.hash = path
  }

  // Render the appropriate page component based on current path
  const renderPage = () => {
    switch (currentPath) {
      case '/':
        return <Dashboard />
      case '/explorer':
        return <ExplorerPage />
      case '/predict':
        return <PredictPage />
      case '/model':
        return <ModelPage />
      case '/docs':
        return <DocsPage />
      default:
        return <Dashboard />
    }
  }

  return (
    <AppProvider>
      <Layout currentPath={currentPath} onNavigate={navigate}>
        <Suspense fallback={
          <div className="flex items-center justify-center min-h-screen">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-yellow-500"></div>
          </div>
        }>
          {renderPage()}
        </Suspense>
      </Layout>
    </AppProvider>
  )
}

export default App
