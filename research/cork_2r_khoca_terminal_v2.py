#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import cork_2r_khoca_terminal as core

extra = json.loads(
    (Path(__file__).with_name("cork_2r_extra_candidates.json")).read_text()
)
for name, record in extra.items():
    core.CANDIDATES[name] = record["pd"]

if __name__ == "__main__":
    core.main()
