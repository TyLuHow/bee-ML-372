"""Simple aggressive Unicode fixer"""
import re
import os

def fix_file(filepath):
    """Remove all non-ASCII characters from print statements."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        fixed_lines = []
        changed = False
        
        for line in lines:
            original = line
            # If line contains print, remove non-ASCII
            if 'print(' in line or 'print (' in line:
                # Keep only ASCII printable characters
                line = ''.join(c if ord(c) < 128 else ' ' for c in line)
            fixed_lines.append(line)
            if line != original:
                changed = True
        
        if changed:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.writelines(fixed_lines)
            return True
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False

files = [
    'src/temporal_analysis.py',
    'src/chemical_space.py',
    'src/toxicophores.py',
    'src/recommendations.py',
    'src/source_comparison.py',
    'src/preprocessing.py',
    'run_all_analyses.py',
    'test_scaffold_split.py',
    'run_comprehensive_analysis.py',
]

for f in files:
    if os.path.exists(f):
        if fix_file(f):
            print(f"Fixed: {f}")
        else:
            print(f"OK: {f}")

