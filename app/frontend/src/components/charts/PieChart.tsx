import React from 'react'
import {
  PieChart as RechartsPieChart,
  Pie,
  Cell,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from 'recharts'
import { PieChartData } from '../../types/explorer'

interface PieChartProps {
  data: PieChartData[]
  colors?: string[]
  title?: string
  height?: number
}

const DEFAULT_COLORS = [
  'hsl(262, 80%, 50%)',   // Purple
  'hsl(31, 100%, 50%)',   // Orange
  'hsl(188, 100%, 35%)',  // Teal
  'hsl(324, 100%, 50%)',  // Magenta
  'hsl(45, 100%, 75%)',   // Light honey
  'hsl(142, 76%, 36%)',   // Green
  'hsl(0, 84%, 60%)',     // Red
  'hsl(217, 91%, 60%)',   // Blue
]

export const PieChart: React.FC<PieChartProps> = ({
  data,
  colors = DEFAULT_COLORS,
  title,
  height = 300,
}) => {
  return (
    <div className="w-full">
      {title && (
        <h4 className="text-lg font-semibold text-text-primary mb-4">{title}</h4>
      )}
      <ResponsiveContainer width="100%" height={height}>
        <RechartsPieChart>
          <Pie
            data={data}
            dataKey="value"
            nameKey="name"
            cx="50%"
            cy="50%"
            outerRadius={80}
            label={({ name, percent }) => `${name} ${(percent * 100).toFixed(0)}%`}
            labelLine
          >
            {data.map((entry, index) => (
              <Cell
                key={`cell-${index}`}
                fill={entry.color || colors[index % colors.length]}
              />
            ))}
          </Pie>
          <Tooltip
            contentStyle={{
              backgroundColor: 'white',
              border: '1px solid #e5e7eb',
              borderRadius: '8px',
              padding: '8px',
            }}
          />
          <Legend
            verticalAlign="bottom"
            height={36}
          />
        </RechartsPieChart>
      </ResponsiveContainer>
    </div>
  )
}
