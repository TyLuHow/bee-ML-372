# ApisTox Testing Guide

Quick reference for testing the ML backend and frontend integration.

## Local Testing Commands

### Start Backend
```bash
cd backend
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn api.main:app --reload --port 8000
```

### Start Frontend
```bash
# In separate terminal, from project root
npm install
npm run dev
```

## API Testing

### Health Check
```bash
curl http://localhost:8000/health
```

**Expected Response:**
```json
{
  "status": "healthy",
  "model_loaded": true,
  "model_name": "Random Forest"
}
```

### Test Endpoint (Sample Prediction)
```bash
curl http://localhost:8000/test
```

This tests Imidacloprid (known neonicotinoid insecticide).

### Manual Prediction Request

**Toxic Compound (Imidacloprid):**
```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "smiles": "C1=CN=C(N1)NC(=O)NCCl",
    "category": "Insecticide",
    "mw": 255.66,
    "logP": 0.57,
    "exposure": "Contact, Oral"
  }'
```

**Safe Compound (Glyphosate):**
```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Glyphosate",
    "category": "Herbicide",
    "mw": 169.07,
    "logP": -3.4,
    "exposure": "Contact"
  }'
```

**With SMILES (Atrazine):**
```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Atrazine",
    "smiles": "CCNc1nc(NC(C)C)nc(Cl)n1",
    "category": "Herbicide",
    "mw": 215.68,
    "logP": 2.61,
    "exposure": "Oral"
  }'
```

**Invalid SMILES (Error Handling):**
```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Invalid",
    "smiles": "NOT_A_VALID_SMILES",
    "category": "Insecticide",
    "mw": 200,
    "logP": 2.0,
    "exposure": "Contact"
  }'
```

## Test Cases

### 1. Known Toxic Compounds

| Compound | Category | Expected Result |
|----------|----------|----------------|
| Imidacloprid | Insecticide | Toxic (>80% confidence) |
| Clothianidin | Insecticide | Toxic (>75% confidence) |
| Thiamethoxam | Insecticide | Toxic (>75% confidence) |
| Fipronil | Insecticide | Toxic (>70% confidence) |

### 2. Known Safe Compounds

| Compound | Category | Expected Result |
|----------|----------|----------------|
| Glyphosate | Herbicide | Safe (>75% confidence) |
| 2,4-D | Herbicide | Safe (>70% confidence) |
| Mancozeb | Fungicide | Safe (>70% confidence) |

### 3. Edge Cases

| Test Case | Expected Behavior |
|-----------|------------------|
| Missing SMILES | Uses default molecular descriptors, lower confidence |
| Invalid SMILES | Returns "Uncertain" with error message |
| Missing required fields | Returns 422 validation error |
| Very high LogP (>8) | May indicate unusual compound |
| Very low LogP (<-5) | May indicate ionic compound |

## Frontend Testing

### Browser Console Checks

1. Open DevTools (F12) → Console
2. Submit a compound analysis
3. Look for:
   ```
   Fetching: POST http://localhost:8000/predict
   Response: {toxicity: "Toxic", confidence: 89.3, ...}
   ```

### Network Tab Checks

1. Open DevTools → Network tab
2. Submit analysis
3. Find `/predict` request
4. Verify:
   - Status: 200
   - Request payload matches API schema
   - Response structure correct
   - No CORS errors

### UI Testing Scenarios

**Scenario 1: Standard Insecticide**
1. Enter compound name: "Test Insecticide"
2. Select category: Insecticide
3. MW: 250, LogP: 3.5
4. Exposure: Contact
5. Submit → Verify result displays

**Scenario 2: Herbicide with SMILES**
1. Add SMILES: "CCNc1nc(NC(C)C)nc(Cl)n1"
2. Category: Herbicide
3. Submit → Verify molecular descriptors computed

**Scenario 3: Backend Unavailable**
1. Stop backend server
2. Submit analysis
3. Verify fallback message appears
4. Verify app still functional

## Production Testing

### Railway Backend

```bash
# Replace with your Railway URL
RAILWAY_URL="https://your-app.railway.app"

# Health check
curl $RAILWAY_URL/health

# Test prediction
curl -X POST $RAILWAY_URL/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Glyphosate",
    "category": "Herbicide",
    "mw": 169.07,
    "logP": -3.4,
    "exposure": "Contact"
  }'
```

### Vercel Frontend

1. Visit Vercel deployment URL
2. Open DevTools
3. Submit analysis
4. Verify API calls to Railway URL
5. Check for CORS errors
6. Verify predictions display correctly

## Performance Testing

### Latency Benchmarks

**Expected response times:**
- `/health`: <50ms
- `/predict` (first request): <1000ms (cold start)
- `/predict` (subsequent): <200ms
- Frontend API call: <500ms total

### Load Testing (Optional)

```bash
# Install Apache Bench
# Ubuntu: apt-get install apache2-utils
# Mac: brew install ab

# Test 100 requests, 10 concurrent
ab -n 100 -c 10 -T 'application/json' \
  -p request.json \
  http://localhost:8000/predict
```

Where `request.json` contains a test compound.

## Debugging

### Backend Not Starting

**Check:**
```bash
cd backend
python3 -c "import fastapi, uvicorn, sklearn, xgboost, rdkit"
```

If imports fail:
```bash
pip install -r requirements.txt
```

### Model Loading Errors

**Check model files:**
```bash
ls -lh backend/models/
# Should show: best_model_random_forest.pkl, best_model_xgboost.pkl, preprocessor.pkl
```

**Test model loading:**
```bash
cd backend
python3 -c "
import pickle
with open('models/best_model_random_forest.pkl', 'rb') as f:
    model = pickle.load(f)
print('Model loaded:', type(model))
"
```

### CORS Errors

**Update backend/api/main.py:**
```python
allow_origins=[
    "http://localhost:3000",
    "http://localhost:5173",
    "https://your-specific-vercel-url.vercel.app",
]
```

### Frontend Not Connecting

**Check environment:**
```bash
# Verify VITE_API_URL is set
cat .env.local

# Should show:
# VITE_API_URL=http://localhost:8000
```

**Test manually:**
```javascript
// In browser console
fetch('http://localhost:8000/health')
  .then(r => r.json())
  .then(console.log)
```

## Common SMILES for Testing

### Insecticides
- Imidacloprid: `C1=CN=C(N1)NC(=O)NCCl`
- Clothianidin: `C1=CN=C(N1)NC(=O)N(CCl)C`
- Thiamethoxam: `CN1COCN(C1=N[N+](=O)[O-])Cc2cnc(s2)Cl`
- Fipronil: `C1=C(C(=C(C(=C1Cl)NC(=O)C(C(F)(F)F)S(=O)C(F)(F)F)Cl)C#N)C(F)(F)F`

### Herbicides
- Glyphosate: `C(C(=O)O)NCP(=O)(O)O`
- Atrazine: `CCNc1nc(NC(C)C)nc(Cl)n1`
- 2,4-D: `Clc1cc(Cl)c(OCC(O)=O)cc1`

### Fungicides
- Mancozeb: `C2H4NS2Mn.xZn`
- Azoxystrobin: `COC=C(C(=O)OC)c1ccccc1-c2cc(C)c(Oc3ccc(cn3)C#N)nc2`

## Validation Checklist

Before declaring complete:

- [ ] Backend starts without errors
- [ ] Health endpoint returns healthy status
- [ ] Test endpoint returns prediction
- [ ] Manual predictions work for toxic compound
- [ ] Manual predictions work for safe compound
- [ ] SMILES validation works
- [ ] Invalid SMILES returns error
- [ ] Frontend connects to backend
- [ ] Frontend displays predictions
- [ ] Fallback works when backend unavailable
- [ ] CORS configured correctly
- [ ] Environment variables set
- [ ] Documentation complete

## Support

If tests fail:
1. Check backend logs for errors
2. Verify all dependencies installed
3. Check model files present
4. Verify CORS configuration
5. Test with `/test` endpoint
6. Review DEPLOYMENT_GUIDE.md

## Interactive Testing

Visit http://localhost:8000/docs for Swagger UI where you can:
- View all endpoints
- Test API interactively
- See request/response schemas
- Download OpenAPI spec
