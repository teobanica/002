"""Reproduce both exact C-sector searches and all certificates from source."""
from pathlib import Path
import subprocess,sys
R=Path(__file__).resolve().parents[1];OLD=R.parent/'n8_aaa_complete/code'
def run(args):
 print('RUN',' '.join(map(str,args)),flush=True)
 subprocess.run(list(map(str,args)),check=True)
def py(name):run([sys.executable,R/'code'/name])
def compile(name,source):
 run(['g++','-O3','-std=c++17','-ftrapv','-I',OLD,source,'-o',R/'checks'/name])
def engine(name,output,*args):
 run([sys.executable,OLD/'run_exact.py','--log',R/'checks'/f'{name}.log',R/'outputs'/output,'--',R/'checks'/name,*args])
compile('c_pair_rows',R/'code/c_pair_rows.cpp');engine('c_pair_rows','c_pair_rows.json')
compile('independent_c_rows',R/'code/independent_c_rows.cpp');engine('independent_c_rows','independent_c_rows.json')
py('verify_reduction.py');py('prepare_fixed_pairs.py')
compile('fixed_pair_compatibility',OLD/'generic_pairs.cpp');engine('fixed_pair_compatibility','fixed_pair_compatibility.json',R/'inputs/fixed_pair_families.txt')
py('prepare_cliques.py')
compile('c_pair_cliques',R/'code/c_pair_cliques.cpp');engine('c_pair_cliques','c_square_loci.json',R/'inputs/fixed_pair_cliques.txt')
py('verify_c_squares.py')
compile('canonical_squares',OLD/'canonical_squares.cpp');engine('canonical_squares','c_canonical.json',R/'inputs/canonical_squares.txt')
compile('independent_fixed_squares',R/'code/independent_fixed_squares.cpp');engine('independent_fixed_squares','independent_fixed_squares.json',R/'inputs/fixed_pair_families.txt')
py('verify_final.py')
print('ALL C-SECTOR CHECKS PASSED.',flush=True)
