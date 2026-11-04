"""Rebuild the complete generic AAA calculation from the archived triple input.
Run from any working directory. Requires g++ and SymPy 1.14.0.
Every command is checked before its dependent command starts.
"""
from pathlib import Path
import subprocess,sys,json
R=Path(__file__).resolve().parents[1]
def run(*args):subprocess.run([str(x) for x in args],cwd=R,check=True)
def py(script,*args):run(sys.executable,R/'code'/script,*args)
def enumerate_json(exe,inp,out,n):
    py('run_exact.py','--records',n,'--log',R/'checks'/f'rebuild_{exe}.log',R/'outputs'/out,'--',R/'checks'/exe,R/'inputs'/inp)
for name in ['generic_full','generic_pairs','family_cliques','family_equivalence']:
    run('g++','-O3','-std=c++17','-ftrapv',R/'code'/f'{name}.cpp','-o',R/'checks'/name)
enumerate_json('generic_full','components.txt','generic_all.json',38)
py('prepare_pipeline_inputs.py','maximal')
enumerate_json('generic_pairs','maximal_families.txt','maximal_pairs.json',2)
py('verify_maximal.py')
py('reduce_row_families.py')
enumerate_json('generic_pairs','reduced_families.txt','reduced_pairs.json',38)
py('prepare_pipeline_inputs.py','cliques')
enumerate_json('family_cliques','reduced_cliques.txt','generic_square_loci.json',38)
py('verify_generic_squares.py')
py('certify_families.py')
py('prepare_equivalence.py','--expanded')
enumerate_json('family_equivalence','equivalence_expanded.txt','generic_equivalences_expanded.json',1519)
py('verify_equivalences.py')
print('Complete generic calculation rebuilt and verified.')
