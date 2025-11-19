import React from 'react'
import { Card, CardContent } from '../ui/Card'

interface StatsCardProps {
  title: string
  value: string | number
  subtitle?: string
  icon?: React.ReactNode
  color?: string
  className?: string
}

export const StatsCard: React.FC<StatsCardProps> = ({
  title,
  value,
  subtitle,
  icon,
  color = 'text-bee-yellow',
  className = '',
}) => {
  return (
    <Card className={`hover:shadow-md transition-shadow ${className}`}>
      <CardContent className="p-6">
        <div className="flex items-start justify-between">
          <div className="flex-1">
            <p className="text-sm font-medium text-text-tertiary uppercase tracking-wide mb-2">
              {title}
            </p>
            <p className={`text-3xl font-bold ${color} mb-1`}>
              {value}
            </p>
            {subtitle && (
              <p className="text-sm text-text-secondary">
                {subtitle}
              </p>
            )}
          </div>
          {icon && (
            <div className={`text-4xl ${color} opacity-20`}>
              {icon}
            </div>
          )}
        </div>
      </CardContent>
    </Card>
  )
}
