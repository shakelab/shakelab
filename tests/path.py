try:
    from pathlib import Path
    ROOT = Path(__file__).parent[1]
    DATA = ROOT / 'data'

except:
    import os
    ROOT = os.path.dirname(os.path.abspath(__file__))
    DATA = os.path.join(ROOT, 'data')

print(DATA)


