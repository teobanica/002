"""Rebuild and verify the CCC classification from scratch."""
from pathlib import Path
from collections import Counter
import json,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
start=time.perf_counter()
for script in ('independent_enumeration.py','coverage.py','independent_verification.py','check_special_examples.py'):
    print('\nRUN',script,flush=True)
    subprocess.run([sys.executable,str(ROOT/script)],check=True)
e=json.loads((ROOT/'enumeration.json').read_text())
c=json.loads((ROOT/'coverage.json').read_text())
v=json.loads((ROOT/'independent_verification.json').read_text())
xy=json.loads((ROOT/'special_example_checks.json').read_text())
assert e['counts']['systems']==3780
assert e['counts']['consistent_systems']==1478
assert e['counts']['raw_connected_components']==1876
assert e['distinct_normalized_components']==len(c['certificates'])==838
assert len(c['targets'])==17
assert Counter(len(t['family']['pivot']) for t in c['targets'])=={3:1,2:2,1:3,0:11}
assert v['total_permutations']==514080
assert xy['status']=='passed' and xy['source_cases']==24
assert xy['permutations_checked']==725760
summary={'status':'passed','systems':3780,'independent_solvers':2,
         'consistent_systems':1478,'connected_branches':1876,
         'distinct_normalized_tori':838,'continuous_family_dimensions':[3,2,2,1,1,1],
         'isolated_CCC_classes':11,'equivalence_permutations_checked':514080,
         'special_example_choices_verified':24,'special_example_permutations_checked':725760}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('\nALL CHECKS PASSED',json.dumps(summary),flush=True)
print('Elapsed',round(time.perf_counter()-start,2),'seconds.',flush=True)

