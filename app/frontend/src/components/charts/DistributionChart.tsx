import React from 'react'
import {
  BarChart as RechartsBarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from 'recharts'
import { DistributionData } from '../../types/explorer'

interface DistributionChartProps {
  data: DistributionData[]
  xLabel: string
  yLabel: string
  title?: string
  height?: number
}

export const DistributionChart: React.FC<DistributionChartProps> = ({
  data,
  xLabel,
  yLabel,
  title,
  height = 400,
}) => {
  return (
    <div className="w-full">
      {title && (
        <h4 className="text-lg font-semibold text-text-primary mb-4">{title}</h4>
      )}
      <ResponsiveContainer width="100%" height={height}>
        <RechartsBarChart data={data} margin={{ top: 10, right: 30, left: 0, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          <XAxis
            dataKey="bin_center"
            label={{ value: xLabel, position: 'insideBottom', offset: -10 }}
            tick={{ fontSize: 12 }}
          />
          <YAxis
            label={{ value: yLabel, angle: -90, position: 'insideLeft' }}
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
          <Legend
            wrapperStyle={{ paddingTop: '20px' }}
            iconType="rect"
          />
          <Bar
            dataKey="toxic"
            fill="hsl(0, 84%, 60%)"
            name="Toxic"
            radius={[4, 4, 0, 0]}
          />
          <Bar
            dataKey="non_toxic"
            fill="hsl(142, 76%, 36%)"
            name="Non-Toxic"
            radius={[4, 4, 0, 0]}
          />
        </RechartsBarChart>
      </ResponsiveContainer>
    </div>
  )
}
