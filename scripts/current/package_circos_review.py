from circos_review_common import *
LIMITS=constrain_process()
import sys, hashlib, zipfile, shutil
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib_cache')
sys.path.insert(0,str(OUT/'vendor'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
stems=['A_lovastatin_regions_circos','B_all_genomes_coordinate_circos','C_all_genomes_qc_evidence_circos']
fig,axes=plt.subplots(1,3,figsize=(18,7))
for ax,stem in zip(axes,stems):
    ax.imshow(plt.imread(OUT/'qa'/(stem+'_pdf.png')));ax.axis('off')
fig.subplots_adjust(left=.005,right=.995,top=.96,bottom=.01,wspace=.015)
fig.suptitle('Terreus | Three circular designs for review',fontsize=17,color='#213c50',x=.025,ha='left')
fig.savefig(OUT/'review_contact_sheet.png',dpi=150,facecolor='white');plt.close(fig)
status=readjson(OUT/'data/validation.json');status['visual_review']='Final PDF renders inspected: all three circles, genome key, full linear overview and top/middle/bottom linear crops; no clipping or overlapping labels observed.'
status['report_filter_test']='Case-insensitive search and reset passed in Node VM with a DOM stub; this was not a browser screenshot test.'
writejson(OUT/'data/validation.json',status)
shutil.copy2(BASE/'package_circos_review.py',OUT/'scripts/package_circos_review.py')
files=[p for p in OUT.rglob('*') if p.is_file() and not any(x in p.relative_to(OUT).parts for x in ('vendor','matplotlib_cache','__pycache__')) and p.name!='file_manifest.json']
manifest={str(p.relative_to(OUT)).replace('\\','/'):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)}
writejson(OUT/'file_manifest.json',manifest);files.append(OUT/'file_manifest.json')
archive=OUT.parent/'terreus_circos_review_2026-11-13.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for path in files:z.write(path,str(Path(OUT.name)/path.relative_to(OUT)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print('PACKAGE',archive,archive.stat().st_size,'bytes;',len(files),'files; archive integrity passed')
