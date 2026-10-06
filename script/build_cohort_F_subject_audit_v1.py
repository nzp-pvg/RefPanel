#!/usr/bin/env python3
import csv
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
data_root = repo_root / 'intermediate_data'
p = data_root / 'cohort_F' / 'derived' / 'GSE243474_sample_manifest_v1.tsv'
out = data_root / 'cohort_F' / 'derived' / 'GSE243474_subject_audit_v1.tsv'
rows=list(csv.DictReader(p.open(),delimiter='\t'))
med=[r for r in rows if r['assay_inferred']=='MeDIP']
groups={}
for r in med:
 g=r['characteristics'].split('cell type: ')[1].split(';')[0]; groups.setdefault((g,r['subject_proxy']),[]).append(r)
outrows=[]
for r in med:
 g=r['characteristics'].split('cell type: ')[1].split(';')[0]; key=(g,r['subject_proxy']); rr=sorted(groups[key],key=lambda x:x['geo_accession']); first=rr[0]['geo_accession']==r['geo_accession']
 r=dict(r); r['candidate_subject']=r['subject_proxy']; r['collapse_status']='KEEP_PROXY_BASELINE' if first else 'DUPLICATE_OR_REPEAT_REVIEW'; r['primary_eligible']='PENDING_METADATA_CONFIRMATION'; outrows.append(r)
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=outrows[0].keys(),delimiter='\t');w.writeheader();w.writerows(outrows)
print('MeDIP records',len(med),'proxy units',len(groups),'duplicate/repeat review',sum(r['collapse_status']!='KEEP_PROXY_BASELINE' for r in outrows))
