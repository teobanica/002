"""Reproduce A/B certificate checks, optionally all exhaustive enumerations.

Run in an extracted copy of the complete archive. --full can take many hours.
The existing positive equivalence certificates are verified; rediscovering them
is unnecessary. Their search implementations are included separately.
"""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import argparse,sys,subprocess,json,os
R=Path(__file__).resolve().parents[1];CODE=R/'code';SHARED=R.parent/'n8_aaa_complete/code'
p=argparse.ArgumentParser();p.add_argument('--full',action='store_true');p.add_argument('--final-only',action='store_true');a=p.parse_args()
assert not(a.full and a.final_only)

def run(name,*args):
 print('RUN',name,*map(str,args),flush=True)
 subprocess.run([sys.executable,str(CODE/name),*map(str,args)],check=True,cwd=R.parent)
def output(cmd,path):
 print('RUN',*map(str,cmd),flush=True)
 proc=subprocess.run(list(map(str,cmd)),stdout=subprocess.PIPE,check=True,cwd=R.parent)
 json.loads(proc.stdout);pending=path.with_suffix('.pending');pending.write_bytes(proc.stdout);os.replace(pending,path)
def compile_one(name):
 subprocess.run(['g++','-std=c++17','-O3','-ftrapv','-I',str(SHARED),str(CODE/(name+'.cpp')),'-o',str(R/'checks'/name)],check=True)
def seed_text(path,cs):
 with path.open('w') as f:
  f.write(str(len(cs))+'\n')
  for c in cs:
   f.write(f"{c.get('id',c['class_id'])} {c.get('branch',0)} {c['dim']}\n")
   for row in c['rows']:f.write(' '.join(str(v) for e in row for v in e)+'\n')
def finish_stage(W,fast=False):
 mod=json.loads((W/'outputs/strata.json').read_text())['modulus']
 output([R/'checks/full_pairs',W/'inputs/pairs_strata.txt',mod],W/'outputs/pairs_strata.json')
 run('reduce_relations.py',W);run('prepare_cliques.py',W);run('run_strata.py',W)
 run('extract_squares_fast.py' if fast else 'extract_squares.py',W)

if a.full:
 for name in ['full_rows','full_rows_general','full_pairs','full_cliques','exceptional_frames']:compile_one(name)
 run('prepare_circles.py')
 output([R/'checks/full_rows',R/'inputs/circle_components.txt'],R/'outputs/full_rows_all.json')
 run('reduce_rows.py');run('prepare_strata.py');finish_stage(R)
 cs=json.loads((R.parent/'proofs/AAA_classes11_19/inputs/components.json').read_text())
 inp=R/'inputs/components_11_19.txt';seed_text(inp,cs)
 output([R/'checks/full_rows',inp],R/'outputs/full_rows_11_19.json')
 run('stratify_11_19.py');run('prepare_11_19.py');finish_stage(R/'stage11_19')
 run('run_uniform_rows.py');run('stratify_uniform.py');run('reduce_uniform_circles.py')
 # A positive alternate-frame certificate suffices. Regenerate the complete
 # frame-search input from the actual strata before using it.
 W=R/'stage_uniform';ss=json.loads((W/'outputs/strata_equivalent.json').read_text())['strata']
 cs=json.loads((R.parent/'proofs/consolidation/outputs/canonical.json').read_text())['triples']['AAA']['classes']
 targets=[c for c in cs if c['id'] in (11,19,21,22)]
 # Include all positive-dimensional strata and all non-aliased points.
 chosen=[s for s in ss if s['base']>0 or 'alias' not in s];inp=W/'inputs/exceptional_frames.txt'
 with inp.open('w') as f:
  f.write(str(len(targets))+'\n')
  for c in targets:f.write(str(c['id'])+' '+' '.join(map(str,c['tensor']))+'\n')
  f.write(str(len(chosen))+'\n')
  for s in chosen:
   mod=lcm(30,*(F(*e[0]).denominator for row in s['rows'] for e in row))
   f.write(f"{s['id']} {s['base']} {mod}\n")
   for row in s['rows']:f.write(' '.join(str(v) for e in row for v in [int(F(*e[0])*mod)]+e[1:])+'\n')
 output([R/'checks/exceptional_frames',inp],W/'outputs/exceptional_frames.json')
 run('prune_uniform.py');run('prepare_stratified_stage.py',W,'strata_pruned.json');finish_stage(W,True)
 run('prepare_bbb.py');run('run_bbb_rows.py');run('stratify_bbb.py')
 W=R/'stage_bbb';run('prepare_stratified_stage.py',W,'strata_bbb.json');finish_stage(W)
 doc=json.loads((W/'inputs/reduced_components.json').read_text());mod=doc['modulus'];inp=W/'inputs/raw_pairs.txt'
 with inp.open('w') as f:
  f.write(str(len(doc['components']))+'\n')
  for c in doc['components']:
   data=json.loads((W/f"outputs/full_rows_{c['id']:02d}.json").read_text())[0]
   f.write(f"{c['id']} {c['dim']} {len(data['families'])}\n")
   for row in data['families']:
    f.write(str(row['free'])+'\n'+' '.join(str(v) for e in row['row'] for v in e)+'\n'+str(len(row['eq']))+'\n')
    for e in row['eq']:f.write(' '.join(map(str,e))+'\n')
 output([R/'checks/full_pairs',inp,mod],W/'outputs/raw_pairs.json')
if not a.final_only:
 run('verify_rows.py');run('verify_uniform_rows.py');run('verify_aliases.py');run('prune_uniform.py')
 run('verify_smith.py',R/'stage_uniform');run('verify_bbb.py')
run('final_verification.py')
print('All requested A/B reproduction checks passed.',flush=True)
