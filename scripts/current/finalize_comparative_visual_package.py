from pathlib import Path
import json,hashlib,shutil,subprocess,sys,zipfile,datetime
BASE=Path(__file__).resolve().parent;P=BASE/'comparative_visualization/2026-11-11';PRO=P/'provenance';PRO.mkdir(exist_ok=True);snap=PRO/'scripts';snap.mkdir(exist_ok=True);OLD=BASE/'taxonomy_reference_build/2026-11-10/analysis'
scripts=['fetch_lov_reference.py','prepare_lov_genomic_screen.py','run_lov_genomic_screen.py','run_lov_genomic_screen_capped.sh','summarize_lov_genomic_screen.py','rescue_lov_proteins.py','rescue_lov_proteins_capped.sh','build_lovastatin_visual_data.py','plot_lovastatin_comparisons.py','annotate_terreus_trees_v2.py','terreus_tree_helpers.py','audit_genome_type_anchors.py','prepare_core66_marker_combinations.py','infer_core66_combinations.py','infer_core66_combinations_capped.sh','prepare_reference_pair_trees.py','infer_reference_pair_trees.py','infer_reference_pair_trees_capped.sh','render_core66_and_pair_trees.py','compare_terreus_tree_topologies.py','compare_all_core66_combinations.py','assess_core66_bootstrap_RF.py','build_comparative_visual_report.py','validate_comparative_visual_package.py','finalize_comparative_visual_package.py']
scripts.append('refresh_terreus_annotation_legend.py')
copies=[]
for name in scripts:
    f=BASE/name;assert f.exists();shutil.copy2(f,snap/name);copies.append({'script':name,'original_path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(PRO/'script_manifest.json').write_text(json.dumps(copies,indent=2))
inputs=P/'sources/input_trees';inputs.mkdir(exist_ok=True);aln=P/'sources/reference_alignments';aln.mkdir(exist_ok=True);im=[]
for name in ['BenA','CaM','RPB2','combined_all_strains']:
    for ext in ['treefile','iqtree','contree','ufboot']:
        src=OLD/'trees'/(name+'.'+ext);dst=inputs/src.name;shutil.copy2(src,dst);im.append({'source':str(src),'copy':dst.relative_to(P).as_posix(),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
for g in ['BenA','CaM','RPB2']:
    src=OLD/'tree_inputs'/(g+'_all_strains.fasta');dst=aln/src.name;shutil.copy2(src,dst);im.append({'source':str(src),'copy':dst.relative_to(P).as_posix(),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
src=Path(r'D:\Illumina\David_Eisele\terreus_genotyping\results\0_4_mafft_iqtree\Sep28\concatenated.fa.treefile');dst=inputs/'core_genome_original.treefile';shutil.copy2(src,dst);im.append({'source':str(src),'copy':dst.relative_to(P).as_posix(),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
for name in ['tree_reference_metadata.json','genome_sources.json','local_metadata.json','sample_taxonomic_review.json']:
    src=OLD/name;dst=P/'sources'/name;shutil.copy2(src,dst);im.append({'source':str(src),'copy':dst.relative_to(P).as_posix(),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
(PRO/'input_manifest.json').write_text(json.dumps(im,indent=2))
(PRO/'execution_notes.txt').write_text('New analysis package: 11 November 2026.\nSource project D:\\Illumina\\David_Eisele\\terreus_genotyping was read only. The November 10 reference package was reused without edits.\n\nLinux sequence stages were launched serially with:\nwsl -d Ubuntu -- systemd-run --user --unit=<unique> --wait --pipe -p MemoryMax=2G -p MemoryHigh=1G -p MemorySwapMax=0 -p CPUQuota=100% -p Nice=15 -p IOWeight=10 /bin/bash <stage-script>\nEach stage also sets ulimit -v 1572864 and one numerical worker/thread.\nServices: terreus-lov-screen-20261111, terreus-core66-combinations-20261111, terreus-lov-protein-rescue-20261111, terreus-reference-pairs-20261111.\nPanaroo was not modified or signalled. The cap is configuration evidence, not a claimed measured peak. Windows plotting and metadata steps were separate.\n\nBLAST+ 2.12.0+; IQ-TREE 2.4.0; scikit-learn 1.8.0. Existing marker alignment source: MAFFT 7.526, 10 November validated package.\nExact sequence commands: lovastatin/blastn_execution.json, lovastatin/tblastn/execution.json, core66_combinations/execution.json, reference_pair_combinations/execution.json.\n\nType-clade rings use external curated ex-type references. Reduced trees have separate projected-reference and documented-genome-type evidence rings. Three core genomes are accepted culture type anchors; six clinical type-material claims are withheld for identity review.\n\nSource snapshots retain their original absolute/workspace-relative paths. They are exact reviewable code, not an automatically portable workflow. Whole source genomes are not copied; extracted genomic context windows and input manifests are included.\n',encoding='utf-8')
(P/'START_HERE.txt').write_text('Open terreus_comparative_figures.html for the figure explorer, all seven marker comparisons, type-clade keys and the 77-genome lovastatin atlas.\n\nCore66: independently inferred and reference-pruned versions, all exact 66 tips. All seven marker choices plus original core-genome topology compared by unrooted RF, NMDS and metric MDS.\n\nLovastatin: 77 assembly labels / 72 distinct canonical sequence sets. Linear and circular views with core-only and up-to-20-kb surroundings. Repeated aliases are retained; unknown is not absence.\n\nFigures include editable SVG, PDF and PNG. An all-genome 77-page PDF atlas is in figures/lovastatin/. The original source project and Panaroo were not modified.\n',encoding='utf-8')
subprocess.run([sys.executable,str(BASE/'validate_comparative_visual_package.py')],check=True)
files=[]
for f in sorted(P.rglob('*')):
    if f.is_file() and f.name!='file_manifest.json':files.append({'path':f.relative_to(P).as_posix(),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(P/'file_manifest.json').write_text(json.dumps(files,indent=2))
archive=P.parent/'terreus_comparative_figures_2026-11-11.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=4) as z:
    for f in sorted(P.rglob('*')):
        if f.is_file():z.write(f,'terreus_comparative_figures_2026-11-11/'+f.relative_to(P).as_posix())
with zipfile.ZipFile(archive) as z:assert z.testzip() is None;assert len(z.namelist())==len(files)+1
assert all(hashlib.sha256((P/f['path']).read_bytes()).hexdigest()==f['sha256'] for f in files)
print('FINAL PACKAGE',archive,'MiB',round(archive.stat().st_size/1024**2,1),'FILES',len(files)+1)
