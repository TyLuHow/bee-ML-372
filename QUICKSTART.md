# ApisTox Quick Start Guide

Get the ML-powered bee toxicity prediction API running in 5 minutes.

## Prerequisites

- Python 3.9+ (`python3 --version`)
- Node.js 18+ (`node --version`)
- pip (`python3 -m pip --version`)

## 1. Start Backend (Terminal 1)

```bash
cd backend
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn api.main:app --reload --port 8000
```

**Wait for:**
```
INFO:     Loaded Random Forest model from ...
INFO:     Application startup complete.
INFO:     Uvicorn running on http://0.0.0.0:8000
```

## 2. Verify Backend

Open new terminal:
```bash
curl http://localhost:8000/health
```

**Expected:**
```json
{
  "status": "healthy",
  "model_loaded": true,
  "model_name": "Random Forest"
}
```

## 3. Test Prediction

```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "category": "Insecticide",
    "mw": 255.66,
    "logP": 0.57,
    "exposure": "Contact"
  }'
```

**Expected:**
```json
{
  "toxicity": "Toxic",
  "confidence": 85.3,
  "explanation": "Analysis indicates...",
  "recommendation": "Recommend avoiding..."
}
```

## 4. Start Frontend (Terminal 2)

```bash
# From project root
npm install
npm run dev
```

**Wait for:**
```
  VITE v5.x.x  ready in xxx ms

  ➜  Local:   http://localhost:5173/
```

## 5. Test Frontend

1. Open http://localhost:5173
2. Enter compound data:
   - Name: "Test Compound"
   - Category: Insecticide
   - MW: 250
   - LogP: 3.0
   - Exposure: Contact
3. Click "Analyze Compound"
4. Verify prediction displays

## 6. View API Docs

Open http://localhost:8000/docs for interactive Swagger UI

## Common Issues

### Backend won't start
```bash
# Check Python version
python3 --version  # Should be 3.9+

# Install dependencies
cd backend
pip install -r requirements.txt
```

### Frontend can't connect
```bash
# Verify .env.local
cat .env.local
# Should show: VITE_API_URL=http://localhost:8000

# Check backend is running
curl http://localhost:8000/health
```

### Port already in use
```bash
# Use different port
uvicorn api.main:app --reload --port 8001

# Update .env.local
VITE_API_URL=http://localhost:8001
```

## Next Steps

- Read [README.md](README.md) for full documentation
- See [TESTING.md](TESTING.md) for test cases
- Check [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) for production deployment
- View [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md) for technical details

## Quick Test Commands

```bash
# Health check
curl http://localhost:8000/health

# Test endpoint (Imidacloprid)
curl http://localhost:8000/test

# Safe compound (Glyphosate)
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"name":"Glyphosate","category":"Herbicide","mw":169.07,"logP":-3.4,"exposure":"Contact"}'

# With SMILES (Atrazine)
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"name":"Atrazine","smiles":"CCNc1nc(NC(C)C)nc(Cl)n1","category":"Herbicide","mw":215.68,"logP":2.61,"exposure":"Oral"}'
```

## File Structure

```
apistox-pro/
├── backend/              # Python ML API
│   ├── api/             # FastAPI endpoints
│   ├── models/          # Trained ML models (3.6 MB)
│   └── requirements.txt
├── services/            # Frontend API client
├── .env.local           # Environment config
└── README.md            # Full documentation
```

## Support

If you encounter issues:
1. Check backend logs in terminal
2. Check browser console for frontend errors
3. Verify environment variables in `.env.local`
4. Review [TESTING.md](TESTING.md) for debugging steps

## Success!

You should now have:
- ✅ Backend running on http://localhost:8000
- ✅ Frontend running on http://localhost:5173
- ✅ API docs at http://localhost:8000/docs
- ✅ ML predictions working end-to-end
