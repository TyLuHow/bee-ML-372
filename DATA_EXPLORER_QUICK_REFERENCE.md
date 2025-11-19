# Data Explorer - Quick Reference Guide

## Tab Overview

### 1. Overview Tab
**Purpose**: Dataset statistics at a glance
**What you'll see**:
- 6 metric cards showing key statistics
- 2 pie charts for chemical types and data sources
- Dataset description

**Key Metrics**:
- Total compounds: 1,035
- Temporal range: 1960-2023
- Toxic percentage: ~42%
- Data sources: Multiple scientific databases

---

### 2. Molecular Diversity Tab
**Purpose**: Explore molecular descriptor distributions
**What you'll see**:
- 6 histograms comparing toxic vs non-toxic
- Statistics tables below each chart

**Descriptors**:
1. Molecular Weight (MW)
2. LogP (lipophilicity)
3. TPSA (polar surface area)
4. H-bond donors
5. H-bond acceptors
6. Rotatable bonds

**Color coding**:
- Red bars = Toxic compounds
- Green bars = Non-toxic compounds

---

### 3. Chemical Space Tab
**Purpose**: Visualize high-dimensional chemical space
**What you'll see**:
- Large scatter plot (PCA or t-SNE)
- Method toggle buttons
- Color-by dropdown selector

**Controls**:
- **Method**: PCA (global) or t-SNE (local clusters)
- **Color by**: Toxicity / Year / Chemical Type / LogP quartiles

**Interpretation**:
- Red points = Toxic
- Green points = Non-toxic
- Clusters = Similar compounds
- Outliers = Unique structures

---

### 4. Temporal Trends Tab
**Purpose**: Analyze toxicity over time
**What you'll see**:
- Line chart: toxicity rate by decade
- Bar chart: sample counts
- Statistical test results box

**Statistical Test**:
- Mann-Kendall trend test
- Shows if toxicity is increasing/decreasing over time
- Significance markers: *, **, *** (p < 0.05, 0.01, 0.001)

---

### 5. Toxicophores Tab
**Purpose**: Identify toxic structural patterns
**What you'll see**:
- Horizontal bar chart: enrichment ratios
- Scatter plot: prevalence vs toxicity
- Top 10 toxicophores table

**Key Concept**:
- **Toxicophore**: Chemical substructure associated with toxicity
- **Enrichment Ratio**: How much more common in toxic vs non-toxic
- **High enrichment + high significance = strong alert**

**Use Case**: Avoid these patterns when designing new chemicals

---

### 6. Correlations Tab
**Purpose**: Explore feature relationships
**What you'll see**:
- Large correlation heatmap (15×15)
- Top positive correlations list
- Top negative correlations list
- Full table of top 20

**Color Scale**:
- Blue = Negative correlation
- White = No correlation
- Red = Positive correlation

**Interpretation**:
- |r| > 0.7 = Strong correlation (features redundant)
- |r| < 0.2 = Weak correlation (features independent)

**Interactive**: Click cells to see details

---

### 7. Properties Tab
**Purpose**: Explore pairwise property relationships
**What you'll see**:
- 3 scatter plots side-by-side
  1. LogP vs Molecular Weight
  2. TPSA vs H-bond Donors
  3. Rotatable Bonds vs Fraction CSP3

**Color coding**:
- Red = Toxic compounds
- Green = Non-toxic compounds

**Look for**:
- Clustering of toxic compounds
- Separation between toxic/non-toxic
- Violations of drug-like rules

---

## Common Actions

### Navigate Between Tabs
- Click tab buttons at the top
- Tab switching is instant (no reload)

### View Details
- **Hover** over any chart element for tooltip
- **Tooltip shows**: Compound name, values, toxicity status

### Interpret Charts
- Red/toxic = Associated with bee toxicity
- Green/non-toxic = Safe for bees
- Larger scatter points = More samples

### Handle Errors
- If data fails to load, click "Retry" button
- Check internet connection
- Verify backend API is running

---

## Pro Tips

1. **Start with Overview** to understand dataset composition
2. **Use Chemical Space** to identify clusters and outliers
3. **Check Correlations** before building models (avoid redundant features)
4. **Toxicophores** are great for early screening
5. **Properties Tab** helps understand drug-likeness
6. **Temporal Trends** reveals potential biases in data collection

---

## Color Legend

### Toxicity Colors
- 🔴 **Red** (hsl(0, 84%, 60%)) = Toxic
- 🟢 **Green** (hsl(142, 76%, 36%)) = Non-toxic / Safe

### Chart Colors
- 🟣 **Purple** (hsl(262, 80%, 50%)) = Primary
- 🟠 **Orange** (hsl(31, 100%, 50%)) = Secondary
- 🔵 **Teal** (hsl(188, 100%, 35%)) = Tertiary
- 🟡 **Yellow** (#FDB813) = Bee accent

### Statistical Significance
- *** = p < 0.001 (very significant)
- ** = p < 0.01 (significant)
- * = p < 0.05 (marginally significant)
- ns = not significant

---

## Keyboard Shortcuts

- **Tab**: Navigate between tabs
- **Enter/Space**: Activate button
- **Arrow keys**: Navigate within charts
- **Hover**: Show tooltips

---

## Troubleshooting

### Charts not loading?
1. Check browser console (F12)
2. Verify API is running
3. Click "Retry" button

### Slow performance?
1. Close other browser tabs
2. Refresh page
3. Check network speed

### Layout broken?
1. Resize browser window
2. Clear browser cache
3. Try different browser

---

## Data Quality Indicators

### Good Quality
- ✅ Balanced toxic/non-toxic ratio
- ✅ Wide temporal range
- ✅ Multiple data sources
- ✅ Diverse chemical types

### Watch Out For
- ⚠️ Very low sample counts in some decades
- ⚠️ Extreme outliers in chemical space
- ⚠️ Very high feature correlations (>0.9)

---

## Integration with Other Pages

### From Predict Page
- After making prediction, explore similar compounds in Chemical Space

### From Model Page
- Check feature importances, then view correlations
- Verify temporal bias in training data

### To Documentation
- Technical details in DATA_EXPLORER_DOCUMENTATION.md
- API details in backend documentation

---

## Export Options (Coming Soon)

Planned features:
- 💾 Export chart data to CSV
- 🖼️ Download charts as PNG/SVG
- 📄 Generate PDF reports

---

## Support

**For technical issues**:
- Check DATA_EXPLORER_DOCUMENTATION.md
- Review browser console errors
- Contact development team

**For scientific questions**:
- Refer to original papers cited in dataset
- Consult toxicology resources
- Review molecular descriptor definitions

---

## Version Information

- **Version**: 1.0.0
- **Release Date**: November 19, 2025
- **Last Updated**: November 19, 2025
- **Status**: Production Ready

---

**Quick Reference Last Updated**: November 19, 2025
