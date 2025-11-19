import React from 'react'

interface ConfusionMatrixProps {
  truePositive: number
  falsePositive: number
  trueNegative: number
  falseNegative: number
}

export const ConfusionMatrix: React.FC<ConfusionMatrixProps> = ({
  truePositive,
  falsePositive,
  trueNegative,
  falseNegative,
}) => {
  const total = truePositive + falsePositive + trueNegative + falseNegative
  const maxValue = Math.max(truePositive, falsePositive, trueNegative, falseNegative)

  const getOpacity = (value: number) => {
    return 0.2 + (value / maxValue) * 0.8
  }

  const getPercentage = (value: number) => {
    return ((value / total) * 100).toFixed(1)
  }

  const Cell: React.FC<{
    value: number
    label: string
    color: string
    description: string
  }> = ({ value, label, color, description }) => (
    <div
      className={`relative p-6 rounded-lg border-2 cursor-pointer transition-all hover:scale-105 ${color}`}
      style={{ opacity: getOpacity(value) }}
      title={description}
    >
      <div className="text-center">
        <div className="text-3xl font-bold text-gray-800">{value}</div>
        <div className="text-sm font-semibold mt-1 text-gray-700">{label}</div>
        <div className="text-xs text-gray-600 mt-1">{getPercentage(value)}%</div>
      </div>
    </div>
  )

  return (
    <div className="w-full max-w-2xl mx-auto">
      <div className="mb-4">
        <h4 className="text-lg font-semibold text-center text-gray-800 mb-2">
          Confusion Matrix
        </h4>
        <p className="text-xs text-center text-gray-600">
          Model performance on test set (n={total} samples)
        </p>
      </div>

      <div className="grid grid-cols-3 gap-2 mb-2">
        <div></div>
        <div className="text-center text-sm font-semibold text-gray-700">
          Predicted Non-Toxic
        </div>
        <div className="text-center text-sm font-semibold text-gray-700">
          Predicted Toxic
        </div>
      </div>

      <div className="grid grid-cols-3 gap-2">
        <div className="flex items-center justify-end pr-4 text-sm font-semibold text-gray-700">
          Actual Non-Toxic
        </div>
        <Cell
          value={trueNegative}
          label="True Negative"
          color="bg-green-100 border-green-300"
          description="Correctly predicted as non-toxic"
        />
        <Cell
          value={falsePositive}
          label="False Positive"
          color="bg-red-100 border-red-300"
          description="Incorrectly predicted as toxic (Type I error)"
        />

        <div className="flex items-center justify-end pr-4 text-sm font-semibold text-gray-700">
          Actual Toxic
        </div>
        <Cell
          value={falseNegative}
          label="False Negative"
          color="bg-orange-100 border-orange-300"
          description="Incorrectly predicted as non-toxic (Type II error)"
        />
        <Cell
          value={truePositive}
          label="True Positive"
          color="bg-green-100 border-green-300"
          description="Correctly predicted as toxic"
        />
      </div>

      <div className="mt-6 grid grid-cols-2 gap-4 text-sm">
        <div className="bg-gray-50 p-3 rounded-lg">
          <div className="font-semibold text-gray-700">Accuracy</div>
          <div className="text-2xl font-bold text-gray-800">
            {(((truePositive + trueNegative) / total) * 100).toFixed(1)}%
          </div>
        </div>
        <div className="bg-gray-50 p-3 rounded-lg">
          <div className="font-semibold text-gray-700">Precision</div>
          <div className="text-2xl font-bold text-gray-800">
            {((truePositive / (truePositive + falsePositive)) * 100).toFixed(1)}%
          </div>
        </div>
        <div className="bg-gray-50 p-3 rounded-lg">
          <div className="font-semibold text-gray-700">Recall (Sensitivity)</div>
          <div className="text-2xl font-bold text-gray-800">
            {((truePositive / (truePositive + falseNegative)) * 100).toFixed(1)}%
          </div>
        </div>
        <div className="bg-gray-50 p-3 rounded-lg">
          <div className="font-semibold text-gray-700">Specificity</div>
          <div className="text-2xl font-bold text-gray-800">
            {((trueNegative / (trueNegative + falsePositive)) * 100).toFixed(1)}%
          </div>
        </div>
      </div>
    </div>
  )
}
