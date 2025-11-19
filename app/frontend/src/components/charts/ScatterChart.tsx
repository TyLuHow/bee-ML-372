import React from 'react'
import {
  ScatterChart as RechartsScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  Cell,
} from 'recharts'
import { ScatterChartData } from '../../types/explorer'

interface ScatterChartProps {
  data: ScatterChartData[]
  xLabel?: string
  yLabel?: string
  title?: string
  colorBy?: 'toxicity' | 'category' | 'none'
  height?: number
  showLegend?: boolean
}

export const ScatterChart: React.FC<ScatterChartProps> = React.memo(({
  data,
  xLabel,
  yLabel,
  title,
  colorBy = 'none',
  height = 400,
  showLegend = true,
}) => {
  const getColor = (point: ScatterChartData) => {
    if (colorBy === 'toxicity') {
      return point.toxicity === 1 ? 'hsl(0, 84%, 60%)' : 'hsl(142, 76%, 36%)'
    }
    if (colorBy === 'category' && point.color) {
      return point.color
    }
    return 'hsl(262, 80%, 50%)'
  }

  const CustomTooltip = ({ active, payload }: any) => {
    if (active && payload && payload.length) {
      const data = payload[0].payload
      return (
        <div className="bg-white border border-border-default rounded-lg p-3 shadow-lg">
          <p className="font-semibold text-text-primary mb-2">
            {data.label || data.compound_name || 'Compound'}
          </p>
          <p className="text-sm text-text-secondary">
            {xLabel || 'X'}: {data.x.toFixed(2)}
          </p>
          <p className="text-sm text-text-secondary">
            {yLabel || 'Y'}: {data.y.toFixed(2)}
          </p>
          {data.toxicity !== undefined && (
            <p className="text-sm mt-1">
              <span className={data.toxicity === 1 ? 'text-accent-toxic' : 'text-accent-safe'}>
                {data.toxicity === 1 ? 'Toxic' : 'Non-Toxic'}
              </span>
            </p>
          )}
          {data.compound_id && (
            <p className="text-xs text-text-tertiary mt-1">
              ID: {data.compound_id}
            </p>
          )}
        </div>
      )
    }
    return null
  }

  return (
    <div className="w-full">
      {title && (
        <h4 className="text-lg font-semibold text-text-primary mb-4">{title}</h4>
      )}
      <ResponsiveContainer width="100%" height={height}>
        <RechartsScatterChart margin={{ top: 10, right: 30, left: 0, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          <XAxis
            type="number"
            dataKey="x"
            label={xLabel ? { value: xLabel, position: 'insideBottom', offset: -10 } : undefined}
            tick={{ fontSize: 12 }}
          />
          <YAxis
            type="number"
            dataKey="y"
            label={yLabel ? { value: yLabel, angle: -90, position: 'insideLeft' } : undefined}
            tick={{ fontSize: 12 }}
          />
          <Tooltip content={<CustomTooltip />} />
          {showLegend && colorBy === 'toxicity' && (
            <Legend
              payload={[
                { value: 'Toxic', type: 'circle', color: 'hsl(0, 84%, 60%)' },
                { value: 'Non-Toxic', type: 'circle', color: 'hsl(142, 76%, 36%)' },
              ]}
            />
          )}
          <Scatter data={data}>
            {data.map((entry, index) => (
              <Cell key={`cell-${index}`} fill={getColor(entry)} />
            ))}
          </Scatter>
        </RechartsScatterChart>
      </ResponsiveContainer>
    </div>
  )
})
