import React from 'react'
import { RadialBarChart, RadialBar, PolarAngleAxis, ResponsiveContainer } from 'recharts'

interface RadialGaugeProps {
  value: number // 0-100
  label?: string
  color?: string
  size?: number
}

export const RadialGauge: React.FC<RadialGaugeProps> = ({
  value,
  label = 'Confidence',
  color,
  size = 200,
}) => {
  // Determine color based on value if not provided
  const gaugeColor = color || (value >= 80 ? '#10b981' : value >= 60 ? '#f59e0b' : '#ef4444')

  const data = [
    {
      name: label,
      value: value,
      fill: gaugeColor,
    },
  ]

  return (
    <div className="flex flex-col items-center">
      <ResponsiveContainer width={size} height={size}>
        <RadialBarChart
          cx="50%"
          cy="50%"
          innerRadius="60%"
          outerRadius="90%"
          data={data}
          startAngle={90}
          endAngle={-270}
        >
          <PolarAngleAxis
            type="number"
            domain={[0, 100]}
            angleAxisId={0}
            tick={false}
          />
          <RadialBar
            background={{ fill: '#e5e7eb' }}
            dataKey="value"
            cornerRadius={10}
            fill={gaugeColor}
          />
          <text
            x="50%"
            y="50%"
            textAnchor="middle"
            dominantBaseline="middle"
            className="text-3xl font-bold fill-current text-gray-800"
          >
            {value.toFixed(1)}%
          </text>
        </RadialBarChart>
      </ResponsiveContainer>
      <p className="text-sm text-gray-600 mt-2">{label}</p>
    </div>
  )
}
