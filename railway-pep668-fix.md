# Railway Deployment Fix: PEP 668 Externally-Managed Environment Error

## Problem Summary

The ApisTox backend deployment on Railway was failing due to a PEP 668 "externally-managed-environment" error when attempting to install Python packages via pip.

### Error Message
```
error: externally-managed-environment

× This environment is externally managed
╰─> This command has been disabled as it tries to modify the immutable
    `/nix/store` filesystem.

    To use Python with Nix and nixpkgs, have a look at the online documentation:
    <https://nixos.org/manual/nixpkgs/stable/#python>.

note: If you believe this is a mistake, please contact your Python installation
or OS distribution provider. You can override this, at the risk of breaking your
Python installation or OS, by passing --break-system-packages.
hint: See PEP 668 for the detailed specification.
```

## Root Cause

**PEP 668** (introduced in Python 3.11+) prevents `pip install` from modifying system-managed Python installations to avoid conflicts between system package managers (like apt, dnf, nix) and pip.

Railway uses **Nixpacks** as its build system, which creates an immutable Nix-based Python environment. The `/nix/store` filesystem is read-only by design, so standard `pip install` commands fail.

## Solution: Two-Step Fix

### Step 1: Fixed Undefined Variable Error (COMPLETED)
**File:** `nixpacks.toml`

Changed from:
```toml
[phases.setup]
nixPkgs = ["python311", "pip"]  # ❌ 'pip' was undefined
```

To:
```toml
[phases.setup]
nixPkgs = ["python311", "python311Packages.pip"]  # ✅ Correct Nix package reference
```

### Step 2: Fixed PEP 668 Error (CURRENT FIX)
**File:** `nixpacks.toml`

Changed from:
```toml
[phases.install]
cmds = ["pip install -r backend/requirements.txt"]  # ❌ Blocked by PEP 668
```

To:
```toml
[phases.install]
cmds = ["pip install --break-system-packages -r backend/requirements.txt"]  # ✅ Bypasses PEP 668
```

## Why `--break-system-packages` is Safe Here

The `--break-system-packages` flag sounds dangerous, but it's actually the **correct and safe solution** for containerized deployments:

1. **Ephemeral Containers**: Railway containers are rebuilt from scratch on every deployment
2. **No User System at Risk**: This isn't a developer's local machine - it's a disposable build environment
3. **Immutable Infrastructure**: The container is destroyed and recreated; there's no persistent "system" to break
4. **Standard Practice**: This is the recommended approach for Nix + Python 3.11+ in container environments
5. **Railway-Specific**: Nixpacks + modern Python requires this flag for pip to function

## Alternative Approaches (Not Used)

### Option 1: Virtual Environment
```toml
[phases.install]
cmds = [
  "python -m venv /opt/venv",
  ". /opt/venv/bin/activate && pip install -r backend/requirements.txt"
]
```
**Why not used:** Adds complexity; Railway containers don't need venv isolation

### Option 2: User Install
```toml
[phases.install]
cmds = ["python -m pip install --user -r backend/requirements.txt"]
```
**Why not used:** Can cause PATH issues; less reliable than `--break-system-packages` in Nix

### Option 3: Nix-Native Package Management
```toml
[phases.setup]
nixPkgs = [
  "python311",
  "python311Packages.fastapi",
  "python311Packages.uvicorn",
  # ... all packages
]
```
**Why not used:** Not all packages available in nixpkgs (e.g., rdkit-pypi); harder to maintain

## Dependencies Being Installed

From `backend/requirements.txt`:

- **FastAPI Stack**: fastapi, uvicorn, pydantic
- **ML Libraries**: scikit-learn, xgboost, imbalanced-learn
- **Chemistry**: rdkit-pypi (molecular computing)
- **Data Science**: pandas, numpy
- **Utilities**: python-multipart

## Verification Steps

After pushing this fix, verify deployment success:

1. **Build Logs**: Check Railway dashboard for successful build
   - Look for: "Successfully installed [package list]"
   - No PEP 668 errors

2. **Deployment Status**: Railway shows "Healthy" status
   - Healthcheck endpoint: `/health` returns 200 OK

3. **Test Health Endpoint**:
   ```bash
   curl https://[your-railway-url]/health
   ```
   Expected response:
   ```json
   {
     "status": "healthy",
     "model_loaded": true
   }
   ```

4. **Test Prediction Endpoint**:
   ```bash
   curl -X POST https://[your-railway-url]/predict \
     -H "Content-Type: application/json" \
     -d '{"smiles": "CCO", "model": "xgboost"}'
   ```
   Should return ML prediction results

## Files Modified

1. **`/mnt/c/Users/Tyler Luby Howard/apis_tox_dataset/nixpacks.toml`**
   - Line 9: Added `--break-system-packages` flag to pip install command

## Commit Message

```
Fix PEP 668 externally-managed-environment error in Railway deployment

- Add --break-system-packages flag to pip install command in nixpacks.toml
- Required for Python 3.11+ in Nix-based Railway containers
- Allows pip to install packages in immutable /nix/store environment
- Safe for ephemeral container deployments

Resolves: Railway build failure blocking ML model backend deployment
```

## References

- [PEP 668 – Marking Python base environments as "externally managed"](https://peps.python.org/pep-0668/)
- [Nixpkgs Python Documentation](https://nixos.org/manual/nixpkgs/stable/#python)
- [Railway Nixpacks Documentation](https://nixpacks.com/docs)

## Timeline

- **2025-12-07**: First fix applied (python311Packages.pip)
- **2025-12-07**: Second fix applied (--break-system-packages flag)
- **Status**: Ready for deployment testing

---

**Next Steps:**
1. Commit changes to git
2. Push to trigger Railway rebuild
3. Monitor build logs for success
4. Test deployed endpoints
5. Update this document with deployment URL and final verification results
