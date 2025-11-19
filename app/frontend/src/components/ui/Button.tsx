import React from 'react'

interface ButtonProps extends React.ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: 'primary' | 'secondary' | 'outline' | 'ghost' | 'danger'
  size?: 'sm' | 'md' | 'lg'
  children: React.ReactNode
}

export const Button: React.FC<ButtonProps> = ({
  variant = 'primary',
  size = 'md',
  className = '',
  children,
  disabled,
  ...props
}) => {
  const baseStyles = 'inline-flex items-center justify-center rounded-lg font-medium transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 disabled:opacity-50 disabled:cursor-not-allowed'

  const variantStyles = {
    primary: 'bg-bee-yellow text-bee-black hover:bg-honey-medium focus:ring-bee-yellow shadow-sm',
    secondary: 'bg-bg-tertiary text-text-primary hover:bg-gray-300 focus:ring-gray-400',
    outline: 'border-2 border-border-default text-text-primary hover:bg-bg-secondary focus:ring-gray-400',
    ghost: 'text-text-primary hover:bg-bg-secondary focus:ring-gray-400',
    danger: 'bg-accent-toxic text-white hover:bg-red-600 focus:ring-red-500 shadow-sm',
  }

  const sizeStyles = {
    sm: 'px-3 py-1.5 text-sm',
    md: 'px-4 py-2 text-base',
    lg: 'px-6 py-3 text-lg',
  }

  return (
    <button
      className={`${baseStyles} ${variantStyles[variant]} ${sizeStyles[size]} ${className}`}
      disabled={disabled}
      {...props}
    >
      {children}
    </button>
  )
}
