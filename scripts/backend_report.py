import json
from pathlib import Path
from geocusp.backends import detect
Path('results').mkdir(exist_ok=True); Path('results/backend_report.json').write_text(json.dumps(detect(),indent=2)); print(json.dumps(detect(),indent=2))
