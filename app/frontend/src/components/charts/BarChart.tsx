import React from 'react'
import {
  BarChart as RechartsBarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ErrorBar,
} from 'recharts'
import { BarChartData } from '../../types/explorer'

interface BarChartProps {
  data: BarChartData[]
  orientation?: 'horizontal' | 'vertical'
  color?: string
  showError?: boolean
  title?: string
  xLabel?: string
  yLabel?: string
  height?: number
}

export const BarChart: React.FC<BarChartProps> = ({
  data,
  orientation = 'vertical',
  color = 'hsl(262, 80%, 50%)',
  showError = false,
  title,
  xLabel,
  yLabel,
  height = 400,
}) => {
  return (
    <div className="w-full">
      {title && (
        <h4 className="text-lg font-semibold text-text-primary mb-4">{title}</h4>
      )}
      <ResponsiveContainer width="100%" height={height}>
        <RechartsBarChart
          data={data}
          layout={orientation === 'horizontal' ? 'vertical' : 'horizontal'}
          margin={{ top: 10, right: 30, left: orientation === 'horizontal' ? 120 : 0, bottom: 20 }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          {orientation === 'vertical' ? (
            <>
              <XAxis
                dataKey="name"
                label={xLabel ? { value: xLabel, position: 'insideBottom', offset: -10 } : undefined}
                tick={{ fontSize: 12 }}
              />
              <YAxis
                label={yLabel ? { value: yLabel, angle: -90, position: 'insideLeft' } : undefined}
                tick={{ fontSize: 12 }}
              />
            </>
          ) : (
            <>
              <XAxis
                type="number"
                label={xLabel ? { value: xLabel, position: 'insideBottom', offset: -10 } : undefined}
                tick={{ fontSize: 12 }}
              />
              <YAxis
                type="category"
                dataKey="name"
                label={yLabel ? { value: yLabel, angle: -90, position: 'insideLeft' } : undefined}
                tick={{ fontSize: 12 }}
                width={100}
              />
            </>
          )}
          <Tooltip
            contentStyle={{
              backgroundColor: 'white',
              border: '1px solid #e5e7eb',
              borderRadius: '8px',
              padding: '8px',
            }}
          />
          <Bar
            dataKey="value"
            fill={color}
            radius={[4, 4, 0, 0]}
          >
            {showError && <ErrorBar dataKey="error" width={4} strokeWidth={2} />}
          </Bar>
        </RechartsBarChart>
      </ResponsiveContainer>
    </div>
  )
}
