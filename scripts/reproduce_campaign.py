#!/usr/bin/env python3
"""Version-2 producer. See docs/FRESH_CAMPAIGN.md for the evidence contract."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from mp5d_campaign.cli import main  # noqa: E402
if __name__ == '__main__':
    raise SystemExit(main())
