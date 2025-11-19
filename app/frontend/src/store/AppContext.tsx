import React, { createContext, useContext, useState, ReactNode } from 'react'

export interface PredictionResult {
  prediction: number
  label_text: string
  probability_toxic: number
  probability_non_toxic: number
  confidence: number
  timestamp: string
  compound_name?: string
}

interface UserPreferences {
  theme: 'light' | 'dark'
  showNotifications: boolean
}

interface AppState {
  predictions: PredictionResult[]
  userPreferences: UserPreferences
  addPrediction: (prediction: PredictionResult) => void
  clearPredictions: () => void
  updatePreferences: (preferences: Partial<UserPreferences>) => void
}

const defaultPreferences: UserPreferences = {
  theme: 'light',
  showNotifications: true,
}

const AppContext = createContext<AppState | undefined>(undefined)

export const AppProvider: React.FC<{ children: ReactNode }> = ({ children }) => {
  const [predictions, setPredictions] = useState<PredictionResult[]>([])
  const [userPreferences, setUserPreferences] = useState<UserPreferences>(defaultPreferences)

  const addPrediction = (prediction: PredictionResult) => {
    setPredictions((prev) => [prediction, ...prev].slice(0, 50)) // Keep last 50 predictions
  }

  const clearPredictions = () => {
    setPredictions([])
  }

  const updatePreferences = (preferences: Partial<UserPreferences>) => {
    setUserPreferences((prev) => ({ ...prev, ...preferences }))
  }

  const value: AppState = {
    predictions,
    userPreferences,
    addPrediction,
    clearPredictions,
    updatePreferences,
  }

  return <AppContext.Provider value={value}>{children}</AppContext.Provider>
}

export const useAppContext = () => {
  const context = useContext(AppContext)
  if (context === undefined) {
    throw new Error('useAppContext must be used within an AppProvider')
  }
  return context
}
