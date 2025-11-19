import React from 'react'

interface BadgeProps {
  children: React.ReactNode
  variant?: 'default' | 'toxic' | 'safe' | 'warning' | 'info'
  className?: string
}

export const Badge: React.FC<BadgeProps> = ({
  children,
  variant = 'default',
  className = '',
}) => {
  const baseStyles = 'inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium'

  const variantStyles = {
    default: 'bg-bg-tertiary text-text-primary',
    toxic: 'bg-red-100 text-red-800',
    safe: 'bg-green-100 text-green-800',
    warning: 'bg-yellow-100 text-yellow-800',
    info: 'bg-blue-100 text-blue-800',
  }

  return (
    <span className={`${baseStyles} ${variantStyles[variant]} ${className}`}>
      {children}
    </span>
  )
}
