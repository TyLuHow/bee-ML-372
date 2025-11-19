import React from 'react'
import { Card, CardHeader, CardTitle, CardContent } from '../components/ui/Card'
import { Badge } from '../components/ui/Badge'

export const DocsPage: React.FC = () => {
  return (
    <div className="space-y-6 max-w-4xl">
      {/* Introduction */}
      <Card>
        <CardHeader>
          <CardTitle>Welcome to ApisTox Documentation</CardTitle>
        </CardHeader>
        <CardContent className="prose prose-sm max-w-none">
          <p className="text-text-secondary">
            ApisTox is an ML-powered application for predicting pesticide toxicity to honey bees.
            This documentation will help you understand how to use the application effectively.
          </p>
        </CardContent>
      </Card>

      {/* Getting Started */}
      <Card>
        <CardHeader>
          <CardTitle>Getting Started</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            <div>
              <h4 className="font-semibold text-text-primary mb-2">1. Navigate to Predict Toxicity</h4>
              <p className="text-sm text-text-secondary">
                Click on "Predict Toxicity" in the sidebar to access the prediction form.
              </p>
            </div>
            <div>
              <h4 className="font-semibold text-text-primary mb-2">2. Enter Compound Data</h4>
              <p className="text-sm text-text-secondary mb-2">
                Fill in the molecular descriptors for your pesticide compound. Required fields include:
              </p>
              <ul className="text-sm text-text-secondary space-y-1 ml-4">
                <li>• Molecular weight and formula</li>
                <li>• Lipophilicity (logP)</li>
                <li>• Hydrogen bond donors/acceptors</li>
                <li>• Rotatable bonds</li>
                <li>• Aromatic rings</li>
              </ul>
            </div>
            <div>
              <h4 className="font-semibold text-text-primary mb-2">3. Get Prediction</h4>
              <p className="text-sm text-text-secondary">
                Click "Predict Toxicity" to receive an instant assessment with confidence scores.
              </p>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Understanding Results */}
      <Card>
        <CardHeader>
          <CardTitle>Understanding Results</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-3">
            <div className="p-4 bg-bg-secondary rounded-lg">
              <div className="flex items-center gap-2 mb-2">
                <Badge variant="toxic">Toxic</Badge>
                <span className="text-sm font-medium">High Risk to Bees</span>
              </div>
              <p className="text-sm text-text-secondary">
                The compound is predicted to be harmful to honey bees. Consider alternative pesticides.
              </p>
            </div>
            <div className="p-4 bg-bg-secondary rounded-lg">
              <div className="flex items-center gap-2 mb-2">
                <Badge variant="safe">Non-Toxic</Badge>
                <span className="text-sm font-medium">Low Risk to Bees</span>
              </div>
              <p className="text-sm text-text-secondary">
                The compound is predicted to be safe for honey bees under normal usage conditions.
              </p>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Model Information */}
      <Card>
        <CardHeader>
          <CardTitle>Model Specifications</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 gap-4 text-sm">
            <div>
              <p className="text-text-tertiary">Algorithm</p>
              <p className="font-medium text-text-primary">XGBoost Classifier</p>
            </div>
            <div>
              <p className="text-text-tertiary">Accuracy</p>
              <p className="font-medium text-accent-safe">83.6%</p>
            </div>
            <div>
              <p className="text-text-tertiary">Training Dataset</p>
              <p className="font-medium text-text-primary">1,035 compounds</p>
            </div>
            <div>
              <p className="text-text-tertiary">Interpretability</p>
              <p className="font-medium text-text-primary">SHAP values</p>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* API Reference */}
      <Card>
        <CardHeader>
          <CardTitle>API Reference</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            <div>
              <h4 className="font-semibold text-text-primary mb-2">Prediction Endpoint</h4>
              <div className="bg-bee-black text-white p-4 rounded-lg font-mono text-sm">
                <div className="mb-2">
                  <Badge variant="info">POST</Badge> /api/predict
                </div>
                <pre className="text-xs overflow-x-auto">
{`{
  "features": {
    "molecular_weight": 180.2,
    "logP": 2.5,
    "h_donors": 2,
    "h_acceptors": 3,
    ...
  }
}`}
                </pre>
              </div>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
