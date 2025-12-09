"""Vercel serverless function wrapper for FastAPI backend"""
import sys
import os

# Add backend to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'backend'))

from api.main import app

# Export the FastAPI app for Vercel
handler = app
