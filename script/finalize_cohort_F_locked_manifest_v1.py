import csv
from pathlib import Path
repo_root = Path(__file__).resolve().parents[2]
root = repo_root / 'intermediate_data' / 'cohort_F'
inp=root/'derived/GSE243474_subject_audit_v1.tsv'; out=root/'derived/GSE243474_locked_manifest_v1.tsv'
rows=list(csv.DictReader(inp.open(),delimiter='\t'))
for r in rows:
    r['adjudication']='KEEP_BASELINE' if r['collapse_status']=='KEEP_PROXY_BASELINE' else 'EXCLUDE_REPLICATE_B_SUFFIX'
    r['primary_eligible']='YES' if r['collapse_status']=='KEEP_PROXY_BASELINE' else 'NO'
    r['adjudication_basis']='Title proxy: b suffix denotes paired replicate; retain lexicographically first accession' if r['collapse_status']!='KEEP_PROXY_BASELINE' else 'Title proxy baseline'
fields=list(rows[0])
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
print('records',len(rows),'eligible',sum(r['primary_eligible']=='YES' for r in rows),'excluded',sum(r['primary_eligible']=='NO' for r in rows))
