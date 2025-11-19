"""
Fix Unicode Characters in Python Files
=======================================
Replaces Unicode symbols with ASCII equivalents for Windows console compatibility.
"""

import os
import re

# Mapping of Unicode symbols to ASCII
REPLACEMENTS = {
    '✓': 'OK',
    '✗': 'X',
    '⚠️': 'WARNING:',
    '⚠': 'WARNING:',
    '📊': '[DATA]',
    '📂': '[FILE]',
    '🔬': '[ANALYSIS]',
    '⚡': '[STAT]',
    '🔄': '[PROCESS]',
    '✅': '[OK]',
    '❌': '[FAIL]',
    '⏳': '[PENDING]',
    '🎯': '[TARGET]',
    '📈': '[CHART]',
    '🧬': '[DNA]',
    '≥': '>=',
    '≤': '<=',
    '→': '->',
    '←': '<-',
    '±': '+/-',
    '×': 'x',
    '÷': '/',
    '\u2265': '>=',
    '\u2264': '<=',
    '\u2713': 'OK',
    '\u2717': 'X',
}

def fix_unicode_in_file(filepath):
    """Replace Unicode characters in a file."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original_content = content
        
        # Replace each Unicode character
        for unicode_char, ascii_rep in REPLACEMENTS.items():
            content = content.replace(unicode_char, ascii_rep)
        
        # Also replace Unicode escape sequences
        import re
        # Replace \uXXXX and \UXXXXXXXX patterns
        content = re.sub(r'\\u[0-9a-fA-F]{4}', lambda m: REPLACEMENTS.get(bytes.fromhex(m.group()[2:]).decode('utf-16-be'), m.group()), content, flags=0)
        content = re.sub(r'\\U[0-9a-fA-F]{8}', lambda m: '[EMOJI]', content, flags=0)
        
        # Replace literal emoji/special characters with regex
        content = re.sub(r'[^\x00-\x7F]+', lambda m: REPLACEMENTS.get(m.group(), '[?]'), content)
        
        # Only write if changes were made
        if content != original_content:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            return True
        return False
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    print("="*60)
    print("FIXING UNICODE CHARACTERS")
    print("="*60)
    
    # Files to fix
    files_to_fix = [
        'src/temporal_analysis.py',
        'src/chemical_space.py',
        'run_comprehensive_analysis.py',
        'src/toxicophores.py',
        'src/recommendations.py',
        'src/source_comparison.py',
        'test_scaffold_split.py',
    ]
    
    fixed_count = 0
    for filepath in files_to_fix:
        if os.path.exists(filepath):
            print(f"Processing: {filepath}...")
            if fix_unicode_in_file(filepath):
                print(f"  OK Fixed")
                fixed_count += 1
            else:
                print(f"  - No changes needed")
        else:
            print(f"  WARNING: File not found: {filepath}")
    
    print("\n" + "="*60)
    print(f"OK Fixed {fixed_count} files")
    print("="*60)

if __name__ == "__main__":
    main()

