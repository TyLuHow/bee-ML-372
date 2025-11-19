#!/usr/bin/env python3
"""
Caching Utilities for Explorer Backend
========================================

Provides caching decorators and utilities for expensive operations
like PCA, t-SNE, and correlation matrix computations.

Author: IME 372 Project Team
Date: November 2025
"""

from functools import wraps
from typing import Dict, Any, Optional, Callable
import hashlib
import json
import time
import os


class SimpleCache:
    """
    Simple in-memory cache with TTL support.

    Stores cached results with timestamps and evicts expired entries.
    """

    def __init__(self, ttl_seconds: int = 3600):
        """
        Initialize cache.

        Args:
            ttl_seconds: Time-to-live for cached entries (default: 1 hour)
        """
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._ttl = ttl_seconds

    def get(self, key: str) -> Optional[Any]:
        """
        Retrieve value from cache if not expired.

        Args:
            key: Cache key

        Returns:
            Cached value or None if not found/expired
        """
        if key not in self._cache:
            return None

        entry = self._cache[key]

        # Check if expired
        if time.time() - entry['timestamp'] > self._ttl:
            del self._cache[key]
            return None

        return entry['value']

    def set(self, key: str, value: Any) -> None:
        """
        Store value in cache with current timestamp.

        Args:
            key: Cache key
            value: Value to cache
        """
        self._cache[key] = {
            'value': value,
            'timestamp': time.time()
        }

    def clear(self) -> None:
        """Clear all cached entries."""
        self._cache.clear()

    def size(self) -> int:
        """Get number of cached entries."""
        return len(self._cache)

    def evict_expired(self) -> int:
        """
        Remove all expired entries.

        Returns:
            Number of entries evicted
        """
        current_time = time.time()
        expired_keys = [
            key for key, entry in self._cache.items()
            if current_time - entry['timestamp'] > self._ttl
        ]

        for key in expired_keys:
            del self._cache[key]

        return len(expired_keys)


# Global cache instance for explorer endpoints
_explorer_cache = SimpleCache(ttl_seconds=3600)  # 1 hour TTL


def cached(cache_key_prefix: str = "", ttl_seconds: Optional[int] = None):
    """
    Decorator for caching function results.

    Args:
        cache_key_prefix: Prefix for cache key (function name used if empty)
        ttl_seconds: Override default TTL (uses cache instance TTL if None)

    Returns:
        Decorated function with caching

    Example:
        @cached(cache_key_prefix="pca_results")
        def compute_pca(data):
            # expensive computation
            return result
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Create cache key from function name and arguments
            key_parts = [cache_key_prefix or func.__name__]

            # Add args to key (convert to string)
            if args:
                args_str = str(args)
                key_parts.append(hashlib.md5(args_str.encode()).hexdigest()[:8])

            # Add kwargs to key (sorted for consistency)
            if kwargs:
                kwargs_str = json.dumps(kwargs, sort_keys=True)
                key_parts.append(hashlib.md5(kwargs_str.encode()).hexdigest()[:8])

            cache_key = "_".join(key_parts)

            # Try to get from cache
            cached_value = _explorer_cache.get(cache_key)
            if cached_value is not None:
                return cached_value

            # Compute and cache result
            result = func(*args, **kwargs)
            _explorer_cache.set(cache_key, result)

            return result

        return wrapper
    return decorator


def cache_to_file(file_path: str, force_refresh: bool = False):
    """
    Decorator for caching function results to a JSON file.

    Args:
        file_path: Path to cache file
        force_refresh: If True, ignore cached file and recompute

    Returns:
        Decorated function with file-based caching

    Example:
        @cache_to_file("outputs/cache/chemical_space.json")
        def compute_chemical_space():
            # expensive computation
            return result
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Check if cached file exists and is recent
            if not force_refresh and os.path.exists(file_path):
                try:
                    # Check file age
                    file_age = time.time() - os.path.getmtime(file_path)
                    if file_age < 3600:  # 1 hour
                        with open(file_path, 'r') as f:
                            return json.load(f)
                except Exception as e:
                    print(f"Warning: Error reading cache file {file_path}: {e}")

            # Compute result
            result = func(*args, **kwargs)

            # Save to file
            try:
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                with open(file_path, 'w') as f:
                    json.dump(result, f, indent=2)
            except Exception as e:
                print(f"Warning: Could not cache to file {file_path}: {e}")

            return result

        return wrapper
    return decorator


def get_cache_stats() -> Dict[str, Any]:
    """
    Get cache statistics.

    Returns:
        Dictionary with cache metrics
    """
    evicted = _explorer_cache.evict_expired()

    return {
        'size': _explorer_cache.size(),
        'ttl_seconds': _explorer_cache._ttl,
        'evicted_entries': evicted
    }


def clear_cache() -> None:
    """Clear all cached entries."""
    _explorer_cache.clear()
