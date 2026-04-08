import sys
from pathlib import Path

# Make workflow/scripts importable without installation
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
