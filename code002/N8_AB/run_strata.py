"""Checkpoint exact clique enumerations independently, with bounded parallelism."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,subprocess,hashlib,os
HOME_ROOT=Path(__file__).resolve().parents[1]
R=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HOME_ROOT
strata=json.loads((R/'outputs/strata.json').read_text())['strata']
MOD=json.loads((R/'outputs/strata.json').read_text())['modulus']
exe=HOME_ROOT/'checks/full_cliques';inp=R/'inputs/cliques_strata.txt'
fingerprint=hashlib.sha256(exe.read_bytes()+inp.read_bytes()).hexdigest()
def run(s):
 i=s['id'];out=R/f'outputs/cliques_stratum_{i:02d}.json';sig=R/f'checks/cliques_stratum_{i:02d}.sha256';log=R/f'checks/cliques_stratum_{i:02d}.log'
 if out.exists() and sig.exists() and sig.read_text().strip()==fingerprint:return json.loads(out.read_text())[0]
 with log.open('w') as f:r=subprocess.run([str(exe),str(inp),str(i),str(MOD)],stdout=subprocess.PIPE,stderr=f)
 assert r.returncode==0,(i,r.returncode)
 data=json.loads(r.stdout);assert len(data)==1 and data[0]['component_index']==i
 tmp=out.with_suffix('.pending');tmp.write_bytes(r.stdout);os.replace(tmp,out);sig.write_text(fingerprint+'\n')
 print(i,'finished, nodes',data[0]['nodes'],'loci',data[0]['loci'],flush=True);return data[0]
with ThreadPoolExecutor(max_workers=4) as pool:
 data=list(pool.map(run,strata))
(R/'outputs/square_loci_strata.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
print('All',len(data),'strata finished.',flush=True)
