import React from 'react'
import { Tabs, TabsList, TabsTrigger, TabsContent } from '../components/ui/Tabs'
import { OverviewTab } from '../components/explorer/OverviewTab'
import { MolecularDiversityTab } from '../components/explorer/MolecularDiversityTab'
import { ChemicalSpaceTab } from '../components/explorer/ChemicalSpaceTab'
import { TemporalTrendsTab } from '../components/explorer/TemporalTrendsTab'
import { ToxicophoresTab } from '../components/explorer/ToxicophoresTab'
import { CorrelationsTab } from '../components/explorer/CorrelationsTab'
import { PropertiesTab } from '../components/explorer/PropertiesTab'

export const ExplorerPage: React.FC = () => {
  return (
    <div className="space-y-6">
      {/* Page Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-text-primary mb-2">Data Explorer</h1>
        <p className="text-text-secondary">
          Interactive visualizations and statistical analysis of the ApisTox honey bee toxicity dataset.
          Explore molecular diversity, chemical space, temporal trends, and more.
        </p>
      </div>

      {/* Tabbed Interface */}
      <Tabs defaultValue="overview">
        <TabsList>
          <TabsTrigger value="overview">Overview</TabsTrigger>
          <TabsTrigger value="molecular-diversity">Molecular Diversity</TabsTrigger>
          <TabsTrigger value="chemical-space">Chemical Space</TabsTrigger>
          <TabsTrigger value="temporal-trends">Temporal Trends</TabsTrigger>
          <TabsTrigger value="toxicophores">Toxicophores</TabsTrigger>
          <TabsTrigger value="correlations">Correlations</TabsTrigger>
          <TabsTrigger value="properties">Properties</TabsTrigger>
        </TabsList>

        <TabsContent value="overview">
          <OverviewTab />
        </TabsContent>

        <TabsContent value="molecular-diversity">
          <MolecularDiversityTab />
        </TabsContent>

        <TabsContent value="chemical-space">
          <ChemicalSpaceTab />
        </TabsContent>

        <TabsContent value="temporal-trends">
          <TemporalTrendsTab />
        </TabsContent>

        <TabsContent value="toxicophores">
          <ToxicophoresTab />
        </TabsContent>

        <TabsContent value="correlations">
          <CorrelationsTab />
        </TabsContent>

        <TabsContent value="properties">
          <PropertiesTab />
        </TabsContent>
      </Tabs>
    </div>
  )
}
