import React from 'react'
import {
  LineChart as RechartsLineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  Area,
} from 'recharts'
import { LineChartData } from '../../types/explorer'

interface LineChartProps {
  data: LineChartData[]
  xLabel?: string
  yLabel?: string
  title?: string
  showConfidence?: boolean
  confidenceData?: { x: number | string; lower: number; upper: number }[]
  height?: number
  color?: string
}

export const LineChart: React.FC<LineChartProps> = ({
  data,
  xLabel,
  yLabel,
  title,
  showConfidence = false,
  confidenceData,
  height = 400,
  color = 'hsl(262, 80%, 50%)',
}) => {
  return (
    <div className="w-full">
      {title && (
        <h4 className="text-lg font-semibold text-text-primary mb-4">{title}</h4>
      )}
      <ResponsiveContainer width="100%" height={height}>
        <RechartsLineChart data={data} margin={{ top: 10, right: 30, left: 0, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          <XAxis
            dataKey="x"
            label={xLabel ? { value: xLabel, position: 'insideBottom', offset: -10 } : undefined}
            tick={{ fontSize: 12 }}
          />
          <YAxis
            label={yLabel ? { value: yLabel, angle: -90, position: 'insideLeft' } : undefined}
            tick={{ fontSize: 12 }}
          />
          <Tooltip
            contentStyle={{
              backgroundColor: 'white',
              border: '1px solid #e5e7eb',
              borderRadius: '8px',
              padding: '8px',
            }}
          />
          {data.some(d => d.label) && <Legend />}
          {showConfidence && confidenceData && (
            <Area
              type="monotone"
              dataKey="upper"
              stroke="none"
              fill={color}
              fillOpacity={0.1}
            />
          )}
          <Line
            type="monotone"
            dataKey="y"
            stroke={color}
            strokeWidth={2}
            dot={{ fill: color, r: 4 }}
            activeDot={{ r: 6 }}
            name={data[0]?.label}
          />
        </RechartsLineChart>
      </ResponsiveContainer>
    </div>
  )
}
