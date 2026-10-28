"""Compatibility entry point; implementation lives in validation/."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from validation.validate_current_model import main


if __name__ == '__main__':
    main()
