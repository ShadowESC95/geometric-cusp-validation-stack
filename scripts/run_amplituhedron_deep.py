#!/usr/bin/env python
from __future__ import annotations
import argparse, json
from geocusp.amplituhedron_sim import run_study

p=argparse.ArgumentParser()
p.add_argument('--samples-k2', type=int, default=20000)
p.add_argument('--samples-k1', type=int, default=12000)
p.add_argument('--seed', type=int, default=20260918)
a=p.parse_args()
r=run_study(a.samples_k2,a.samples_k1,a.seed)
print(json.dumps(r.report, indent=2))
