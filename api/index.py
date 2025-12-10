"""
Vercel serverless function for ApisTox ML prediction API
Optimized for Vercel's 50MB limit and 10s timeout
"""
import sys
import os

# Add backend to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'backend'))

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from api.main import app as fastapi_app

# Export for Vercel
app = fastapi_app
