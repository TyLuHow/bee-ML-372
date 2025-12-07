# Quick Verification Guide - Railway Backend Fix

## What Was Fixed

The Railway backend was failing to build due to incorrect Nixpacks configuration. The fix has been pushed and Railway should be rebuilding now.

## Immediate Next Steps (5-10 minutes)

### Step 1: Check Railway Deployment (2 minutes)

1. Open Railway dashboard: https://railway.app
2. Find project: `apis_tox_dataset` or `bee-ML-372`
3. Look for latest deployment
4. Check build logs for success

**Success indicators:**
```
✓ Installing dependencies
✓ Loaded Random Forest model
✓ Application startup complete
```

### Step 2: Get Your Railway URL (1 minute)

In Railway dashboard:
- Go to Settings → Domains
- Copy the URL (looks like: `https://apis-tox-dataset-production.up.railway.app`)

### Step 3: Test Backend (2 minutes)

Replace `<YOUR_RAILWAY_URL>` below with your actual URL:

```bash
# Test 1: Health check
curl https://<YOUR_RAILWAY_URL>/health

# Expected: {"status":"healthy","model_loaded":true}

# Test 2: Prediction
curl -X POST https://<YOUR_RAILWAY_URL>/predict \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Imidacloprid",
    "smiles": "C1=CN=C(N1)NC(=O)NCCl",
    "category": "Insecticide",
    "mw": 255.66,
    "logP": 0.57,
    "exposure": "Contact"
  }'

# Expected: Real ML prediction (not fallback)
```

### Step 4: Update Frontend Environment Variable (3 minutes)

**In Vercel Dashboard:**
1. Go to your project (`apistox-pro`)
2. Settings → Environment Variables
3. Add or update:
   - Key: `VITE_API_URL`
   - Value: `https://<YOUR_RAILWAY_URL>` (no trailing slash!)
   - Environments: All (Production, Preview, Development)
4. Click "Redeploy" to apply changes

### Step 5: Test End-to-End (2 minutes)

1. Open your Vercel deployment: `https://apistox-pro.vercel.app`
2. Submit a compound analysis (try Imidacloprid as test)
3. Check prediction results

**Success = NO "fallback prediction" message**

## Quick Troubleshooting

| Issue | Solution |
|-------|----------|
| Railway build failing | Check build logs, verify nixpacks.toml syntax |
| Health check returns 404 | Backend not started, check deployment logs |
| CORS error in browser | Update backend CORS to include Vercel domain |
| Still showing fallback | Verify VITE_API_URL set in Vercel, redeploy frontend |

## What Changed

### Files Modified
- `nixpacks.toml` - Fixed Python/pip package references
- `railway.json` - Cleaned up configuration
- `backend-connectivity-fix.md` - Full documentation (read this for details)

### The Fix
```toml
# BEFORE (broken):
nixPkgs = ['python311', 'pip']  # ❌ 'pip' undefined

# AFTER (fixed):
nixPkgs = ["python311", "python311Packages.pip"]  # ✅ Correct
```

## Need More Help?

Read the full documentation: `/backend-connectivity-fix.md`

That file includes:
- Complete root cause analysis
- Detailed verification steps
- Troubleshooting guide
- All test cases
- Architecture diagrams
