# 🚀 Deployment Guide: Railway + Vercel

This guide splits the application into two parts for optimal hosting:
1.  **Backend (FastAPI + ML Model)**: Hosted on **Railway** (supports large files/memory).
2.  **Frontend (React)**: Hosted on **Vercel** (fast, global CDN).

---

## Part 1: Deploy Backend to Railway

1.  **Sign Up/Login**: Go to [railway.app](https://railway.app/) and log in with GitHub.
2.  **New Project**:
    *   Click **"New Project"** -> **"Deploy from GitHub repo"**.
    *   Select your repository (`bee-ML-372` or similar).
    *   Click **"Deploy Now"**.
3.  **Configure Build**:
    *   Railway will auto-detect Python.
    *   It will use the `Procfile` in the root: `web: python -m uvicorn app.backend.main:app --host 0.0.0.0 --port $PORT`.
4.  **Wait for Build**:
    *   The build might take 2-3 minutes as it installs libraries (numpy, pandas, etc.).
5.  **Generate Domain**:
    *   Go to the **Settings** tab of your Railway service.
    *   Under **"Networking"**, click **"Generate Domain"**.
    *   Copy this URL (e.g., `https://apis-tox-dataset-production.up.railway.app`).
    *   **Test it**: Open that URL + `/docs` in your browser (e.g., `https://.../docs`). You should see the Swagger UI.

---

## Part 2: Deploy Frontend to Vercel

1.  **Sign Up/Login**: Go to [vercel.com](https://vercel.com/) and log in with GitHub.
2.  **Add New Project**:
    *   Click **"Add New..."** -> **"Project"**.
    *   Import your repository.
3.  **Configure Project**:
    *   **Framework Preset**: Select `Vite`.
    *   **Root Directory**: Click "Edit" and select `app/frontend`. **(Important!)**
4.  **Environment Variables**:
    *   Expand **"Environment Variables"**.
    *   Add Key: `VITE_API_URL`
    *   Add Value: Your Railway URL from Part 1 (e.g., `https://apis-tox-dataset-production.up.railway.app`).
        *   *Note: Do not add a trailing slash `/`.*
5.  **Deploy**:
    *   Click **"Deploy"**.
    *   Vercel will build your React site.
6.  **Final Test**:
    *   Click the domain provided by Vercel.
    *   Try running a prediction. If it works, your full stack app is live! 🚀

---

## Troubleshooting

**Backend (Railway)**
*   **Build Fails?**: Check the "Build Logs". If it complains about missing files, ensure you pushed everything to git.
*   **App Crashes?**: Check "Deploy Logs". If it says "ModuleNotFoundError", check `requirements.txt`.

**Frontend (Vercel)**
*   **404 on Page Load?**: Ensure "Root Directory" was set to `app/frontend`.
*   **Prediction Error?**:
    *   Open Browser Console (F12).
    *   Check if the network request is going to `localhost` or your Railway URL.
    *   If `localhost`, you forgot the `VITE_API_URL` variable. Go to Vercel Settings -> Environment Variables, add it, and **Redeploy**.
