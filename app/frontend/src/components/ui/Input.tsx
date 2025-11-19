import React from 'react'

interface InputProps extends React.InputHTMLAttributes<HTMLInputElement> {
  label?: string
  error?: string
  helperText?: string
}

export const Input: React.FC<InputProps> = ({
  label,
  error,
  helperText,
  className = '',
  ...props
}) => {
  return (
    <div className="w-full">
      {label && (
        <label className="block text-sm font-medium text-text-primary mb-1.5">
          {label}
        </label>
      )}
      <input
        className={`w-full px-3 py-2 border rounded-lg text-text-primary bg-white
          focus:outline-none focus:ring-2 focus:ring-bee-yellow focus:border-transparent
          disabled:bg-bg-secondary disabled:cursor-not-allowed
          ${error ? 'border-accent-toxic' : 'border-border-default'}
          ${className}`}
        {...props}
      />
      {error && (
        <p className="mt-1 text-sm text-accent-toxic">{error}</p>
      )}
      {helperText && !error && (
        <p className="mt-1 text-sm text-text-tertiary">{helperText}</p>
      )}
    </div>
  )
}
