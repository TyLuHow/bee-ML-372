import React, { useState } from 'react'

interface HeatmapChartProps {
  matrix: number[][]
  features: string[]
  title?: string
  colorScale?: 'diverging' | 'sequential'
  onCellClick?: (row: number, col: number, value: number) => void
}

export const HeatmapChart: React.FC<HeatmapChartProps> = React.memo(({
  matrix,
  features,
  title,
  colorScale = 'diverging',
  onCellClick,
}) => {
  const [hoveredCell, setHoveredCell] = useState<{ row: number; col: number } | null>(null)

  const getColor = (value: number): string => {
    if (colorScale === 'diverging') {
      // Blue (negative) -> White (zero) -> Red (positive)
      if (value < 0) {
        const intensity = Math.abs(value) * 100
        return `hsl(217, 91%, ${100 - intensity * 0.4}%)`
      } else {
        const intensity = value * 100
        return `hsl(0, 84%, ${100 - intensity * 0.4}%)`
      }
    } else {
      // Sequential: white to purple
      const intensity = value * 100
      return `hsl(262, 80%, ${100 - intensity * 0.5}%)`
    }
  }

  const cellSize = Math.min(40, 600 / features.length)

  return (
    <div className="w-full">
      {title && (
        <h4 className="text-lg font-semibold text-text-primary mb-4">{title}</h4>
      )}
      <div className="overflow-x-auto">
        <div className="inline-block min-w-full">
          {/* Feature labels - top */}
          <div className="flex" style={{ marginLeft: `${cellSize * 3}px` }}>
            {features.map((feature, index) => (
              <div
                key={`header-${index}`}
                className="text-xs text-text-secondary transform -rotate-45 origin-bottom-left"
                style={{
                  width: `${cellSize}px`,
                  height: `${cellSize * 2}px`,
                  marginLeft: index === 0 ? '0' : `-${cellSize * 0.3}px`,
                }}
              >
                <span className="inline-block" style={{ width: `${cellSize * 2}px` }}>
                  {feature}
                </span>
              </div>
            ))}
          </div>

          {/* Heatmap grid */}
          <div>
            {matrix.map((row, rowIndex) => (
              <div key={`row-${rowIndex}`} className="flex items-center">
                {/* Row label */}
                <div
                  className="text-xs text-text-secondary text-right pr-2"
                  style={{ width: `${cellSize * 3}px` }}
                >
                  {features[rowIndex]}
                </div>
                {/* Cells */}
                {row.map((value, colIndex) => {
                  const isHovered = hoveredCell?.row === rowIndex && hoveredCell?.col === colIndex
                  return (
                    <div
                      key={`cell-${rowIndex}-${colIndex}`}
                      className={`border border-gray-200 cursor-pointer transition-all ${
                        isHovered ? 'ring-2 ring-bee-yellow z-10' : ''
                      }`}
                      style={{
                        width: `${cellSize}px`,
                        height: `${cellSize}px`,
                        backgroundColor: getColor(value),
                      }}
                      onMouseEnter={() => setHoveredCell({ row: rowIndex, col: colIndex })}
                      onMouseLeave={() => setHoveredCell(null)}
                      onClick={() => onCellClick?.(rowIndex, colIndex, value)}
                      title={`${features[rowIndex]} vs ${features[colIndex]}: ${value.toFixed(3)}`}
                    >
                      {cellSize > 30 && (
                        <div className="flex items-center justify-center h-full text-xs font-medium">
                          {value.toFixed(2)}
                        </div>
                      )}
                    </div>
                  )
                })}
              </div>
            ))}
          </div>

          {/* Color scale legend */}
          <div className="mt-6 flex items-center justify-center gap-4">
            <span className="text-sm text-text-secondary">
              {colorScale === 'diverging' ? '-1.0' : '0.0'}
            </span>
            <div className="flex h-6 rounded overflow-hidden" style={{ width: '200px' }}>
              {Array.from({ length: 20 }).map((_, i) => {
                const value = colorScale === 'diverging' ? -1 + (i / 19) * 2 : i / 19
                return (
                  <div
                    key={i}
                    style={{
                      width: '10px',
                      backgroundColor: getColor(value),
                    }}
                  />
                )
              })}
            </div>
            <span className="text-sm text-text-secondary">
              {colorScale === 'diverging' ? '+1.0' : '1.0'}
            </span>
          </div>
        </div>
      </div>
    </div>
  )
})
