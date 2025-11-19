# 🚀 START HERE - ApisTox Quick Setup

**Welcome!** This guide gets you from zero to running in ~20 minutes.

---

## ⚡ Super Quick Start (3 Commands)

```bash
# 1. Install everything
pip install -r requirements.txt

# 2. Run automated setup (trains models + runs analyses)
python QUICK_START.bat

# 3. Start the API
python -m uvicorn app.backend.main:app --reload --port 8000
```

Then open http://localhost:8000/docs in your browser!

---

## 📋 What You Need

- ✅ Python 3.9-3.12 installed
- ✅ 8GB RAM minimum
- ✅ 5GB free disk space
- ✅ Internet connection (for pip install)

---

## 🎯 The 5-Minute Version

### Step 1: Install (2 min)
```bash
cd "C:\Users\Tyler Luby Howard\apis_tox_dataset"
pip install -r requirements.txt
```

### Step 2: Verify (30 sec)
```bash
python verify_setup.py
```

### Step 3: Train Models (5-10 min)
```bash
python src/models.py
```

### Step 4: Run Analyses (10 min)
```bash
python run_all_analyses.py
```

### Step 5: Start API (10 sec)
```bash
python -m uvicorn app.backend.main:app --reload --port 8000
```

✅ **Done!** Visit http://localhost:8000/docs

---

## 🖥️ The Full Version (with Frontend)

### Terminal 1 - Backend
```bash
python -m uvicorn app.backend.main:app --reload --port 8000
```

### Terminal 2 - Frontend
```bash
cd app/frontend
npm install
npm run dev
```

✅ **Done!** Visit http://localhost:5173

---

## 🆘 If Something Goes Wrong

```bash
# Check what's missing
python verify_setup.py

# If RDKit fails
pip install rdkit==2023.9.5

# If model not found
python src/models.py

# If port 8000 busy
netstat -ano | findstr :8000
# Then kill the process or use port 8001
```

---

## 📚 Need More Help?

| Question | See |
|----------|-----|
| Detailed setup instructions | `SETUP_GUIDE.md` |
| Command reference | `COMMANDS.md` |
| What got implemented | `IMPLEMENTATION_COMPLETE.md` |
| API documentation | http://localhost:8000/docs |
| Troubleshooting | `SETUP_GUIDE.md` Section 8 |

---

## ✅ Success Checklist

After setup, you should have:

- [ ] ✅ `outputs/models/best_model.pkl` exists
- [ ] ✅ `outputs/models/preprocessor.pkl` exists
- [ ] ✅ 20+ PNG files in `outputs/figures/`
- [ ] ✅ 3 markdown reports in `outputs/`
- [ ] ✅ API running at http://localhost:8000
- [ ] ✅ API docs at http://localhost:8000/docs
- [ ] ✅ Can make predictions via API

**Verify**: Run `python verify_setup.py`

---

## 🎓 For Your IME 372 Project

You now have:

1. ✅ **Trained ML models** (XGBoost, 83.6% accuracy)
2. ✅ **RESTful API** (10+ endpoints)
3. ✅ **React frontend** (interactive UI)
4. ✅ **20+ visualizations** (publication-quality)
5. ✅ **3 analysis reports** (markdown format)
6. ✅ **Complete documentation** (6 guide files)
7. ✅ **Automated testing** (pytest suite)
8. ✅ **Reproducible pipeline** (all scripts included)

**Everything is ready for demonstration and submission!**

---

## 📊 What Each Command Does

```bash
# Trains the main XGBoost model
python src/models.py
# → Creates: outputs/models/best_model.pkl

# Runs temporal analysis + chemical space visualization
python run_comprehensive_analysis.py
# → Creates: 10+ figures, temporal report

# Identifies toxicophores (toxic structural patterns)
python src/toxicophores.py
# → Creates: toxicophore report + figures

# Finds safer alternative compounds
python src/recommendations.py
# → Creates: alternatives.csv + report

# Compares ECOTOX vs PPDB data sources
python src/source_comparison.py
# → Creates: source comparison figure + JSON

# Tests scaffold-based splitting
python test_scaffold_split.py
# → Creates: scaffold split results JSON

# Runs ALL of the above in sequence
python run_all_analyses.py
# → Creates: Everything!
```

---

## 🚀 Launch Commands

```bash
# Start Backend API
python -m uvicorn app.backend.main:app --reload --port 8000

# Start Frontend (in new terminal)
cd app/frontend && npm run dev

# Run All Tests
pytest tests/ -v

# Check System Status
python verify_setup.py
```

---

## 🎯 Your Next 15 Minutes

1. **Install** (2 min): `pip install -r requirements.txt`
2. **Train** (5 min): `python src/models.py`
3. **Analyze** (5 min): `python run_all_analyses.py`
4. **Start** (1 min): `python -m uvicorn app.backend.main:app --reload --port 8000`
5. **Test** (2 min): Visit http://localhost:8000/docs and try predictions

**That's it! You're done!** 🎉

---

## 💡 Pro Tips

- Use `QUICK_START.bat` to automate steps 1-3
- Keep the API terminal open while testing
- Check `outputs/figures/` for all visualizations
- Read reports in `outputs/*.md`
- Use `COMMANDS.md` as a cheat sheet

---

**Questions?** See `SETUP_GUIDE.md` for detailed instructions!

**Ready to deploy?** The system works on Vercel, Railway, or any Python host!

---

*Last updated: November 8, 2025*




