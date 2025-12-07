<objective>
Restore all missing input fields that the ML model requires for accurate toxicity predictions. The design overhaul removed critical features from the dataset that are necessary for the trained model to work properly.
</objective>

<context>
During the UI/UX redesign, many input fields from the original dataset were removed. The current frontend only collects:
- name
- mw (molecular weight)
- logP (lipophilicity)
- exposure (route)
- category (Insecticide/Herbicide/Fungicide/Other)

However, the training dataset at `/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/outputs/dataset_final.csv` contains these fields:
- name, CID, CAS, SMILES
- source (ECOTOX, PPDB, BPDB)
- year
- toxicity_type (Contact, Oral, Systemic, Other)
- herbicide, fungicide, insecticide, other_agrochemical (binary flags 0/1)
- label (target variable)
- ppdb_level

The model was trained on ALL these features plus 15 computed molecular descriptors from SMILES. Missing these fields means the model cannot make accurate predictions.

Current files:
- @App.tsx - Frontend form (missing fields)
- @types.ts - TypeScript interfaces
- @backend/api/main.py - Backend API schema
- @backend/api/predict.py - Prediction logic
- @services/geminiService.ts - API service
</context>

<requirements>
1. **Add ALL missing input fields to the frontend form:**
   - **SMILES** (text input, optional but highly recommended for accuracy)
   - **CAS number** (text input, optional)
   - **Source** (dropdown: ECOTOX, PPDB, BPDB, Unknown)
   - **Year** (number input, default to current year)
   - **Toxicity Type** (dropdown: Contact, Oral, Systemic, Other) - DIFFERENT from exposure route
   - **CID** (number input, optional - PubChem Compound ID)

2. **Keep existing category selection** but also add binary flags:
   - Auto-set herbicide/fungicide/insecticide/other_agrochemical based on category selection
   - Or allow manual override with checkboxes

3. **Update TypeScript interfaces:**
   - Add all new fields to ChemicalData interface in types.ts
   - Mark optional fields as optional
   - Add proper types (string, number, boolean)

4. **Update backend API schema:**
   - Add all new fields to CompoundInput model
   - Mark optional fields appropriately
   - Ensure preprocessor handles these features correctly

5. **Update prediction logic:**
   - Pass all fields to the model
   - Handle missing optional fields gracefully
   - Use SMILES to compute molecular descriptors if provided
   - If SMILES missing, still make prediction with available features

6. **Maintain beautiful UI/UX:**
   - Organize fields into logical sections
   - Use collapsible "Advanced Options" for optional fields
   - Keep the Modern Botanical Atelier aesthetic
   - Add helpful tooltips explaining each field
   - Provide example values

7. **Add data validation:**
   - SMILES validation (if provided)
   - CAS number format validation
   - Year range validation (1800-2030)
   - Required vs optional field handling

8. **Update documentation:**
   - Document what each field means
   - Explain why SMILES is important
   - Provide examples of complete inputs
</requirements>

<implementation>
**Form Organization:**
```
Section 1: Basic Information
- Compound Name (required)
- CAS Number (optional, with format hint)
- SMILES (optional but recommended, with validation)

Section 2: Data Source & Classification
- Source (ECOTOX/PPDB/BPDB/Unknown)
- Year (number, default 2025)
- Toxicity Type (Contact/Oral/Systemic/Other)
- Exposure Route (existing field)

Section 3: Pesticide Classification
- Category buttons (existing: Insecticide/Herbicide/Fungicide/Other)
- Auto-populate binary flags from category
- Optional: Advanced checkboxes for manual override

Section 4: Molecular Properties
- Molecular Weight (existing)
- LogP (existing)
- CID (optional, in collapsed "Advanced" section)
```

**Example Complete Input:**
```javascript
{
  name: "Imidacloprid",
  cas: "138261-41-3",
  smiles: "C1=CN=C(N1)NC(=O)NCCl",
  cid: 86287519,
  source: "PPDB",
  year: 2020,
  toxicity_type: "Contact",
  exposure: "Contact (Direct Spray)",
  category: "Insecticide",
  herbicide: 0,
  fungicide: 0,
  insecticide: 1,
  other_agrochemical: 0,
  mw: 255.66,
  logP: 0.57
}
```

**Backend Changes:**
- Update CompoundInput to accept all fields
- Modify preprocessing to use all available features
- Add SMILES validation
- Handle missing values appropriately

**Frontend Changes:**
- Expand form with new sections
- Add collapsible "Advanced Options"
- Implement field validation
- Update state management
- Maintain responsive design
</implementation>

<output>
Modify the following files:
- `./App.tsx` - Add all missing input fields with proper UI organization
- `./types.ts` - Update ChemicalData interface with all fields
- `./backend/api/main.py` - Update CompoundInput schema
- `./backend/api/predict.py` - Handle all input features
- `./services/geminiService.ts` - Pass all fields to backend
</output>

<verification>
Before declaring complete:
1. Verify all dataset fields are now input fields
2. Test that SMILES input triggers molecular descriptor computation
3. Confirm backend accepts all new fields
4. Test prediction with complete vs minimal input
5. Verify UI remains beautiful and usable
6. Check that binary flags auto-populate from category
7. Test form validation for each field
8. Ensure missing optional fields don't break predictions
</verification>

<success_criteria>
- All 13+ dataset fields available as inputs
- SMILES field with validation
- CAS number field with format hints
- Source, year, toxicity_type dropdowns
- Binary flags (herbicide/fungicide/insecticide/other_agrochemical)
- Form organized into logical sections
- Advanced fields in collapsible section
- Beautiful UI maintained
- Backend handles all fields correctly
- Predictions work with partial or complete data
- Field validation implemented
- Example values provided in tooltips/placeholders
</success_criteria>
