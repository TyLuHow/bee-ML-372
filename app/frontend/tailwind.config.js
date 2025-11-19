/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // Background colors
        'bg-primary': 'hsl(0, 0%, 100%)',
        'bg-secondary': 'hsl(210, 20%, 98%)',
        'bg-tertiary': 'hsl(210, 20%, 95%)',

        // Text colors
        'text-primary': 'hsl(222, 47%, 11%)',
        'text-secondary': 'hsl(215, 16%, 47%)',
        'text-tertiary': 'hsl(216, 12%, 68%)',

        // Accent colors
        'accent-toxic': 'hsl(0, 84%, 60%)',
        'accent-safe': 'hsl(142, 76%, 36%)',
        'accent-warning': 'hsl(38, 92%, 50%)',
        'accent-info': 'hsl(217, 91%, 60%)',

        // Chart colors
        'chart-primary': 'hsl(262, 80%, 50%)',
        'chart-secondary': 'hsl(31, 100%, 50%)',
        'chart-tertiary': 'hsl(188, 100%, 35%)',
        'chart-quaternary': 'hsl(324, 100%, 50%)',

        // Border color
        'border-default': 'hsl(214, 32%, 91%)',

        // Honey/Bee theme colors
        'honey-light': 'hsl(45, 100%, 75%)',
        'honey-medium': 'hsl(38, 100%, 50%)',
        'honey-dark': 'hsl(35, 85%, 35%)',
        'bee-yellow': '#FDB813',
        'bee-black': '#1A1A1A',
      },
      fontFamily: {
        sans: ['Inter', 'system-ui', '-apple-system', 'BlinkMacSystemFont', '"Segoe UI"', 'Roboto', 'sans-serif'],
        mono: ['"JetBrains Mono"', '"Fira Code"', 'Consolas', 'Monaco', 'monospace'],
        display: ['Inter', 'sans-serif'],
      },
      fontSize: {
        xs: '0.75rem',
        sm: '0.875rem',
        base: '1rem',
        lg: '1.125rem',
        xl: '1.25rem',
        '2xl': '1.5rem',
        '3xl': '1.875rem',
        '4xl': '2.25rem',
        '5xl': '3rem',
      },
      fontWeight: {
        normal: 400,
        medium: 500,
        semibold: 600,
        bold: 700,
      },
      lineHeight: {
        tight: 1.25,
        normal: 1.5,
        relaxed: 1.75,
      },
      borderRadius: {
        sm: '0.25rem',
        DEFAULT: '0.5rem',
        md: '0.5rem',
        lg: '0.75rem',
        xl: '1rem',
        '2xl': '1.5rem',
      },
      boxShadow: {
        sm: '0 1px 2px 0 rgb(0 0 0 / 0.05)',
        DEFAULT: '0 1px 3px 0 rgb(0 0 0 / 0.1), 0 1px 2px -1px rgb(0 0 0 / 0.1)',
        md: '0 4px 6px -1px rgb(0 0 0 / 0.1), 0 2px 4px -2px rgb(0 0 0 / 0.1)',
        lg: '0 10px 15px -3px rgb(0 0 0 / 0.1), 0 4px 6px -4px rgb(0 0 0 / 0.1)',
        xl: '0 20px 25px -5px rgb(0 0 0 / 0.1), 0 8px 10px -6px rgb(0 0 0 / 0.1)',
      },
    },
  },
  plugins: [],
}
