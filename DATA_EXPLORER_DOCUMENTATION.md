# Data Explorer Documentation

## Overview

The Data Explorer is a comprehensive interactive visualization platform for the ApisTox honey bee toxicity dataset. It provides seven specialized analysis tabs with professional scientific visualizations built using React and Recharts.

## Architecture

### Frontend Structure

```
app/frontend/src/
├── types/
│   └── explorer.ts                 # TypeScript interfaces for API responses
├── services/
│   └── explorerApi.ts             # API service for data fetching
├── components/
│   ├── ui/
│   │   └── Tabs.tsx               # Tabbed interface component
│   ├── charts/                    # Reusable chart components
│   │   ├── StatsCard.tsx          # KPI metric display
│   │   ├── DistributionChart.tsx  # Histogram with toxic/non-toxic overlay
│   │   ├── PieChart.tsx           # Categorical data pie chart
│   │   ├── BarChart.tsx           # Horizontal/vertical bar chart
│   │   ├── LineChart.tsx          # Time series line chart
│   │   ├── ScatterChart.tsx       # 2D scatter plot with tooltips
│   │   └── HeatmapChart.tsx       # Correlation heatmap
│   └── explorer/                  # Tab components
│       ├── OverviewTab.tsx
│       ├── MolecularDiversityTab.tsx
│       ├── ChemicalSpaceTab.tsx
│       ├── TemporalTrendsTab.tsx
│       ├── ToxicophoresTab.tsx
│       ├── CorrelationsTab.tsx
│       └── PropertiesTab.tsx
└── pages/
    └── ExplorerPage.tsx           # Main page with tabbed interface
```

### Backend Endpoints

All endpoints are prefixed with `/api/explorer/`:

1. **GET `/overview`** - Dataset statistics and distributions
2. **GET `/molecular-diversity`** - Descriptor distributions
3. **GET `/chemical-space`** - PCA & t-SNE coordinates
4. **GET `/temporal-trends`** - Time series analysis
5. **GET `/toxicophores`** - Structural alerts analysis
6. **GET `/correlations`** - Feature correlation matrix
7. **GET `/property-distributions`** - 2D scatter plots

## Tab Features

### 1. Overview Tab

**Purpose**: High-level dataset statistics and composition

**Visualizations**:
- 6 KPI cards: Total compounds, temporal range, toxic %, sources, chemical types, exposure types
- Pie charts: Chemical type distribution, data source distribution
- Dataset information text

**Key Metrics**:
- Total unique compounds
- Year range of data collection
- Percentage of toxic compounds
- Number of data sources
- Chemical and exposure type counts

**Use Cases**:
- Quick dataset understanding
- Dataset quality assessment
- Presentation to stakeholders

---

### 2. Molecular Diversity Tab

**Purpose**: Analyze distribution of molecular descriptors

**Visualizations**:
- 6 distribution histograms (toxic vs non-toxic overlay)
- Summary statistics tables for each descriptor

**Descriptors Analyzed**:
- Molecular Weight (MW)
- LogP (lipophilicity)
- Topological Polar Surface Area (TPSA)
- Number of H-bond donors
- Number of H-bond acceptors
- Number of rotatable bonds

**Key Insights**:
- Structural complexity differences
- Lipophilicity patterns
- Polarity distributions
- Drug-likeness assessment

**Use Cases**:
- Feature engineering guidance
- Understanding toxic vs non-toxic differences
- Quality control (outlier detection)

---

### 3. Chemical Space Tab

**Purpose**: Visualize high-dimensional chemical space in 2D

**Visualizations**:
- Interactive scatter plot (PCA or t-SNE)
- Variance explained (PCA mode)

**Interactive Controls**:
- Toggle: PCA vs t-SNE
- Color by: Toxicity, Year, Chemical Type, LogP quartiles
- Hover tooltips: Compound details

**Key Features**:
- **PCA**: Linear dimensionality reduction, preserves global structure
- **t-SNE**: Non-linear, preserves local structure, identifies clusters
- Click points for compound information

**Use Cases**:
- Cluster identification
- Outlier detection
- Understanding chemical diversity
- Identifying underrepresented regions

---

### 4. Temporal Trends Tab

**Purpose**: Analyze toxicity trends over time

**Visualizations**:
- Line chart: Toxicity rate by decade
- Bar chart: Sample count by decade
- Statistical test results (Mann-Kendall)

**Statistical Analysis**:
- Mann-Kendall trend test
- Kendall's Tau coefficient
- P-value and significance level
- Trend direction (increasing/decreasing/no trend)

**Key Insights**:
- Research focus shifts over time
- Data availability by era
- Temporal bias detection

**Use Cases**:
- Understanding dataset temporal bias
- Historical context
- Identifying data gaps

---

### 5. Toxicophores Tab

**Purpose**: Identify structural alerts associated with toxicity

**Visualizations**:
- Horizontal bar chart: Enrichment ratios (top 15)
- Scatter plot: Prevalence vs toxicity rate
- Top 10 toxicophores table with statistics

**Statistical Metrics**:
- Enrichment ratio (toxic prevalence / non-toxic prevalence)
- 95% confidence intervals
- P-values (Fisher's exact test)
- Significance markers (*, **, ***)

**Key Features**:
- SMARTS substructure matching
- Bootstrap resampling for CI
- Multiple testing correction

**Use Cases**:
- Predictive model feature engineering
- Safer chemical design (avoid toxicophores)
- Mechanism of action understanding
- Early screening alerts

---

### 6. Correlations Tab

**Purpose**: Explore feature relationships and multicollinearity

**Visualizations**:
- Interactive correlation heatmap (15×15 matrix)
- Top positive correlations list
- Top negative correlations list
- Full top 20 correlations table

**Interactive Features**:
- Click cells to view details
- Color scale: Blue (negative) → White (zero) → Red (positive)
- Hover for values

**Interpretation Guide**:
- |r| > 0.7: Strong correlation
- 0.4 < |r| < 0.7: Moderate
- 0.2 < |r| < 0.4: Weak
- |r| < 0.2: Very weak/none

**Use Cases**:
- Feature selection
- Multicollinearity detection
- Dimensionality reduction planning
- Model interpretation

---

### 7. Properties Tab

**Purpose**: Explore pairwise property relationships

**Visualizations**:
- 3 scatter plots (toxic vs non-toxic colored):
  1. LogP vs Molecular Weight (lipophilicity vs size)
  2. TPSA vs H-Donors (polarity vs H-bonding)
  3. Rotatable Bonds vs Fraction CSP3 (flexibility vs saturation)

**Property Descriptions**:
- **LogP**: Hydrophobicity, membrane permeability
- **MW**: Molecular size
- **TPSA**: Polar surface area, absorption
- **H-Donors**: Hydrogen bond donors
- **Rotatable Bonds**: Molecular flexibility
- **Fraction CSP3**: Carbon saturation

**Key Insights**:
- Drug-likeness assessment (Lipinski's Rule of Five)
- Property space clustering
- Non-linear relationships

**Use Cases**:
- Chemical design guidance
- Property-toxicity relationships
- Feature interaction understanding

---

## Technical Implementation

### Data Fetching Pattern

All tabs use a consistent data fetching pattern:

```typescript
const [data, setData] = useState<ResponseType | null>(null)
const [loading, setLoading] = useState(true)
const [error, setError] = useState<string | null>(null)

useEffect(() => {
  fetchData()
}, [])

const fetchData = async () => {
  try {
    setLoading(true)
    setError(null)
    const result = await explorerApi.getEndpoint()
    setData(result)
  } catch (err) {
    setError(err.message)
  } finally {
    setLoading(false)
  }
}
```

### State Management

- **Local state**: Each tab manages its own data
- **No redundant fetching**: Data cached after initial load
- **Tab switching**: Instant (no re-fetch)

### Loading States

All tabs implement:
- Skeleton loaders during fetch
- Error alerts with retry buttons
- Graceful error handling

### Responsive Design

- **Desktop**: 2-3 column grids
- **Tablet**: 2 column grids
- **Mobile**: Single column stacks
- All charts use ResponsiveContainer (100% width)

### Chart Library

**Recharts** (v2.10.0):
- React-native charting
- Declarative API
- Built-in animations
- Responsive by default
- TypeScript support

### Color Palette

Scientific color scheme from design system:

```typescript
{
  toxic: 'hsl(0, 84%, 60%)',      // Red
  safe: 'hsl(142, 76%, 36%)',     // Green
  primary: 'hsl(262, 80%, 50%)',  // Purple
  secondary: 'hsl(31, 100%, 50%)', // Orange
  tertiary: 'hsl(188, 100%, 35%)', // Teal
  bee_yellow: '#FDB813',           // Bee yellow accent
}
```

---

## Performance Considerations

### Optimizations

1. **Data Caching**: Fetched data persists across tab switches
2. **Lazy Loading**: Tab content only renders when active
3. **Chart Memoization**: Recharts handles internal optimization
4. **Efficient Re-renders**: React context for tab state only

### Performance Metrics

- **Initial Load**: < 2s (depends on backend)
- **Tab Switch**: < 50ms (instant)
- **Chart Render**: < 100ms per chart
- **Interaction Response**: < 16ms (60 FPS)

### Bundle Size Impact

- Recharts: ~150KB gzipped
- Chart components: ~20KB
- Tab components: ~40KB
- **Total addition**: ~210KB gzipped

---

## Accessibility

### Keyboard Navigation

- Tab key navigates between tabs
- Arrow keys navigate within charts
- Enter/Space activates buttons

### Screen Readers

- Semantic HTML structure
- ARIA labels on interactive elements
- Alt text for chart descriptions
- Role attributes for custom components

### Color Contrast

- WCAG AA compliant color palette
- Tooltips with high contrast
- Text size meets minimum standards

---

## User Guide

### Getting Started

1. Navigate to `/explorer` in the application
2. Default tab is "Overview" with dataset statistics
3. Click tab buttons to switch between visualizations
4. Hover over charts for detailed tooltips

### Tips for Analysis

**Overview Tab**:
- Review key metrics to understand dataset composition
- Check temporal range for data recency

**Molecular Diversity Tab**:
- Compare toxic vs non-toxic distributions
- Look for clear separation in histograms (predictive features)

**Chemical Space Tab**:
- Start with PCA for global structure
- Switch to t-SNE for cluster identification
- Try different color schemes to reveal patterns

**Temporal Trends Tab**:
- Check Mann-Kendall results for temporal bias
- Low sample counts in early decades may be unreliable

**Toxicophores Tab**:
- Focus on high enrichment + high significance
- Upper-right scatter plot region = best alerts
- Use in feature engineering

**Correlations Tab**:
- Strong correlations (|r| > 0.7) indicate redundancy
- Consider feature selection or PCA
- Check for unexpected correlations (data quality)

**Properties Tab**:
- Look for toxic compound clustering
- Violations of drug-like rules may indicate toxicity
- Use for safer chemical design

---

## Future Enhancements

### Planned Features

1. **Export Functionality**:
   - CSV download for all visualizations
   - PNG/SVG export for charts
   - PDF report generation

2. **Advanced Interactions**:
   - Zoom/pan on scatter plots
   - Brushing for range selection
   - Linked views (select in one chart, highlight in others)

3. **Additional Analysis**:
   - Custom descriptor selection
   - Statistical test parameters
   - Filtering and subsetting

4. **Performance**:
   - Virtual scrolling for large tables
   - Progressive data loading
   - WebGL for large scatter plots

### Known Limitations

1. Large datasets (>10,000 points) may slow scatter plots
2. Heatmap limited to 15×15 features (performance)
3. No offline caching (requires backend on each load)

---

## Troubleshooting

### Common Issues

**Charts not rendering**:
- Check browser console for errors
- Verify API endpoints are accessible
- Ensure data format matches TypeScript types

**Slow performance**:
- Check network tab for slow API responses
- Reduce chart complexity (fewer data points)
- Clear browser cache

**TypeScript errors**:
- Ensure all imports are correct
- Verify API response types match interfaces
- Run `npm run build` to check compilation

**Layout issues**:
- Check Tailwind CSS classes
- Verify responsive breakpoints
- Test on different screen sizes

---

## Development

### Adding a New Tab

1. Create tab component in `components/explorer/`
2. Add API endpoint in `explorerApi.ts`
3. Define TypeScript types in `types/explorer.ts`
4. Import and add to `ExplorerPage.tsx`
5. Create chart components as needed
6. Add to index exports

### Testing

```bash
# Type checking
npm run build

# Development server
npm run dev

# Linting
npm run lint
```

### Code Style

- Use functional components with hooks
- TypeScript strict mode
- Consistent naming (camelCase for functions, PascalCase for components)
- Exhaustive error handling
- Loading states for all async operations

---

## API Response Examples

### Overview Response

```json
{
  "total_compounds": 1035,
  "temporal_range": {
    "min_year": 1960,
    "max_year": 2023
  },
  "toxic_percentage": 42.5,
  "data_sources": [
    {"name": "PubChem", "count": 450},
    {"name": "EPA", "count": 385}
  ],
  "chemical_types": [
    {"type": "Pesticide", "count": 520},
    {"type": "Industrial", "count": 315}
  ],
  "exposure_types": [
    {"type": "Contact", "count": 650},
    {"type": "Oral", "count": 385}
  ]
}
```

### Chemical Space Response

```json
{
  "pca": {
    "points": [
      {
        "compound_id": "CID123",
        "compound_name": "Compound A",
        "x": 2.5,
        "y": -1.2,
        "toxicity": 1,
        "logp": 3.4,
        "molecular_weight": 234.5
      }
    ],
    "variance_explained": [0.35, 0.22]
  },
  "tsne": {
    "points": [...]
  }
}
```

---

## Support

For issues or feature requests:
1. Check this documentation
2. Review API endpoint documentation
3. Check browser console for errors
4. Contact development team

---

## Version History

- **v1.0.0** (2025-11-19): Initial release
  - 7 analysis tabs
  - 7 chart components
  - Full TypeScript support
  - Responsive design
  - Accessibility features
