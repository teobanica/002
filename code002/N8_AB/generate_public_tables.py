from pathlib import Path
import json
from collections import defaultdict
R=Path(__file__).resolve().parents[1];W=R/'stage_uniform'
cs=json.loads((W/'inputs/components.json').read_text());groups=defaultdict(list)
for c in cs:groups[c['class_id']].append(c)
lines=[r'\begin{longtable}{crrr}',r'\toprule',r'AAA support class & Base dimension & Components & Raw fourth-row records\\',r'\midrule',r'\endfirsthead',r'\toprule',r'AAA support class & Base dimension & Components & Raw fourth-row records\\',r'\midrule',r'\endhead']
total=0
for cid,items in sorted(groups.items()):
 dims={c['dim'] for c in items};assert len(dims)==1
 n=sum(len(json.loads((W/f"outputs/full_rows_{c['component_index']:02d}.json").read_text())[0]['families']) for c in items);total+=n
 lines.append(f"${cid}$ & ${items[0]['dim']}$ & ${len(items)}$ & ${n:,}$"+r'\\')
lines += [r'\midrule',f'Total & & $32$ & ${total:,}$'+r'\\',r'\bottomrule',r'\end{longtable}']
(R.parent/'code/text_uniform_table.tex').write_text('\n'.join(lines)+'\n');assert total==548310
print('Generated the 19-class fourth-row count table.')
