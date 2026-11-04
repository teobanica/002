#!/usr/bin/env python3
"""Prepare the unchanged proof programs and run saved-certificate checks.

Only this launcher is new. Source paths and input bytes are restored from
manifest.json so that the original mathematical programs need no edits.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parent
CHECKS = {
    'n7-ab': 'N7/code/ab/independent_check.py',
    'n7-ccc': 'N7/code/ccc/independent_verification.py',
    'n7-mixed': 'N7/code/mixed/verify_mixed.py',
    'n7-f4': 'N7/code/f4/verify_candidates.py',
    'n7-isolated': 'N7/code/isolated45/verify.py',
    'n7-affine': 'N7/code/affine4/verify_export.py',
    'n7-completion': 'N7/code/completion/check.py',
    'n7-squares': 'N7/code/square/verify_square_families.py',
    'n8-c-reduction': 'N8/n8_C_sector/code/verify_reduction.py',
    'n8-c-squares': 'N8/n8_C_sector/code/verify_c_squares.py',
    'n8-c-final': 'N8/n8_C_sector/code/verify_final.py',
    'n8-aaa-results': 'N8/n8_aaa_complete/code/verify_results.py',
    'n8-aaa-equivalences': 'N8/n8_aaa_complete/code/verify_equivalences.py',
    'n8-class22-display': 'N8/proofs/AAA_class22/code/verify_displayed_data.py',
    'n8-ab-final': 'N8/n8_AB_completion/code/final_verification.py',
    'n8-cycles': 'N8/code/verify_cycle_censuses.py',
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def child(root, name):
    p = PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or '\\' in name:
        raise ValueError('Unsafe package path: ' + name)
    result = root.joinpath(*p.parts)
    if not result.resolve().is_relative_to(root.resolve()):
        raise ValueError('Path escapes work folder: ' + name)
    return result


def verified(data, item, name):
    if len(data) != item['bytes'] or digest(data) != item['sha256']:
        raise ValueError('Package integrity check failed: ' + name)
    return data


def prepare(work):
    raw = (ROOT / 'manifest.json').read_bytes()
    manifest = json.loads(raw)
    stamp = {'manifest_sha256': digest(raw)}
    marker = work / 'code002-work.json'
    if work.exists():
        if not marker.is_file() or json.loads(marker.read_text()) != stamp:
            raise ValueError('Work folder is not from this package. Choose a new --work path; no files were overwritten.')
        for rel, item in manifest['sources'].items():
            verified(child(work, rel).read_bytes(), item, rel)
        return manifest
    # Validate every archive and source before creating an execution folder.
    for rel, item in manifest['sources'].items():
        verified(child(ROOT, item['file']).read_bytes(), item, item['file'])
    archives = {}
    for item in manifest['data'].values():
        archives.setdefault(item['archive'], {})
    for rel, item in manifest['data'].items():
        archives[item['archive']][rel] = item
    for archive, entries in archives.items():
        with zipfile.ZipFile(child(ROOT, archive)) as z:
            if len(z.namelist()) != len(set(z.namelist())) or set(z.namelist()) != set(entries):
                raise ValueError('Unexpected archive members: ' + archive)
            for rel, item in entries.items():
                child(work, rel)
                verified(z.read(rel), item, rel)
    work.mkdir(parents=True, exist_ok=False)
    for rel in manifest['directories']:
        child(work, rel).mkdir(parents=True, exist_ok=True)
    for rel, item in manifest['sources'].items():
        path = child(work, rel)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(child(ROOT, item['file']).read_bytes())
        if path.suffix == '.sh':
            path.chmod(0o755)
    for archive, entries in archives.items():
        with zipfile.ZipFile(child(ROOT, archive)) as z:
            for rel in entries:
                path = child(work, rel)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(z.read(rel))
    marker.write_text(json.dumps(stamp, indent=2) + '\n')
    print('Prepared', work, flush=True)
    return manifest


def invoke(work, rel, args):
    path = child(work, rel)
    if path.suffix == '.py':
        command = [sys.executable, '-u', str(path)]
    elif path.suffix == '.sh':
        command = ['bash', str(path)]
    else:
        raise ValueError('Use a Python or shell entry point; C++ files are compiled by the reproduction scripts.')
    env = dict(os.environ)
    env.pop('PYTHONOPTIMIZE', None)
    env['PYTHON'] = sys.executable
    return subprocess.run(command + args, cwd=path.parent, env=env).returncode


def main():
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('Run without -O or PYTHONOPTIMIZE: the proof programs use assertions.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, default=ROOT.parent / 'code002_work',
                        help='Separate execution folder (default: ../code002_work). Never overwrite an unrelated folder.')
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('prepare', help='Check file hashes and restore paths, without calculations.')
    commands.add_parser('list', help='List certificate checks, without extracting data.')
    check = commands.add_parser('check', help='Run certificate checks, not the full classification searches.')
    check.add_argument('--only', action='append', choices=sorted(CHECKS), help='May be repeated.')
    run = commands.add_parser('run', help='Run a retained program using its original path.')
    run.add_argument('script')
    run.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command == 'list':
        for name, path in CHECKS.items():
            print(f'{name:24s} {path}')
        return
    work = args.work.expanduser().resolve()
    if work == ROOT or work.is_relative_to(ROOT):
        parser.error('Choose a work folder outside code002 to keep source and generated files separate.')
    manifest = prepare(work)
    if args.command == 'prepare':
        print('Source and data are ready. No mathematical computation was run.')
    elif args.command == 'run':
        rel = args.script.removeprefix('code002/')
        if rel not in manifest['sources']:
            parser.error('Unknown original path. See SOURCE_INDEX.md.')
        tail = args.arguments
        if tail[:1] == ['--']:
            tail = tail[1:]
        raise SystemExit(invoke(work, rel, tail))
    else:
        results = []
        for name in args.only or CHECKS:
            print('\nCHECK', name, flush=True)
            start = time.monotonic()
            code = invoke(work, CHECKS[name], [])
            results.append({'check': name, 'program': CHECKS[name],
                            'exit_code': code, 'seconds': round(time.monotonic()-start, 3)})
            (work / 'validation.json').write_text(json.dumps(results, indent=2)+'\n')
            if code:
                raise SystemExit('FAILED: '+name+'; stopped. See '+str(work/'validation.json'))
        print('\nPASS:', len(results), 'certificate checks. Full searches were not rerun.', flush=True)


if __name__ == '__main__':
    main()
