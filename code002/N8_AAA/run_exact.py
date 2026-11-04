"""Run an enumerator, validate its complete JSON, then atomically save it.

No partially streamed JSON file is ever exposed as a completed output.
"""
import argparse,json,os,stat,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('output');p.add_argument('--records',type=int);p.add_argument('--log');p.add_argument('command',nargs=argparse.REMAINDER)
# Options precede the output path, followed by -- and the actual command.
a=p.parse_args();cmd=a.command
if cmd and cmd[0]=='--':cmd=cmd[1:]
if not cmd:raise SystemExit('missing command')
exe=Path(cmd[0])
if exe.is_file():exe.chmod(exe.stat().st_mode|stat.S_IXUSR)
r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
if a.log:
    q=Path(a.log);q.parent.mkdir(parents=True,exist_ok=True);tmp=q.with_name(q.name+'.pending');tmp.write_bytes(r.stderr);os.replace(tmp,q)
if r.returncode:
    sys.stderr.buffer.write(r.stderr);raise SystemExit(r.returncode)
d=json.loads(r.stdout)
if a.records is not None:assert isinstance(d,list) and len(d)==a.records,(a.output,len(d),a.records)
q=Path(a.output);q.parent.mkdir(parents=True,exist_ok=True);tmp=q.with_name(q.name+'.pending');tmp.write_bytes(r.stdout);os.replace(tmp,q)
print('Saved complete JSON:',q)
