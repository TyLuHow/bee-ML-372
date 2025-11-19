export const colors = {
  background: {
    primary: 'hsl(0, 0%, 100%)',      // Pure white
    secondary: 'hsl(210, 20%, 98%)',  // Off-white
    tertiary: 'hsl(210, 20%, 95%)',   // Light gray
  },
  text: {
    primary: 'hsl(222, 47%, 11%)',    // Near black
    secondary: 'hsl(215, 16%, 47%)',  // Medium gray
    tertiary: 'hsl(216, 12%, 68%)',   // Light gray
  },
  accent: {
    toxic: 'hsl(0, 84%, 60%)',        // Vibrant red
    safe: 'hsl(142, 76%, 36%)',       // Forest green
    warning: 'hsl(38, 92%, 50%)',     // Amber
    info: 'hsl(217, 91%, 60%)',       // Blue
  },
  chart: {
    primary: 'hsl(262, 80%, 50%)',    // Purple
    secondary: 'hsl(31, 100%, 50%)',  // Orange
    tertiary: 'hsl(188, 100%, 35%)',  // Teal
    quaternary: 'hsl(324, 100%, 50%)', // Magenta
  },
  border: 'hsl(214, 32%, 91%)',
  honey: {
    light: 'hsl(45, 100%, 75%)',      // Light honey
    medium: 'hsl(38, 100%, 50%)',     // Medium honey/amber
    dark: 'hsl(35, 85%, 35%)',        // Dark honey
  },
  bee: {
    yellow: '#FDB813',                 // Classic bee yellow
    black: '#1A1A1A',                  // Bee black
  }
} as const

export type ColorScheme = typeof colors
