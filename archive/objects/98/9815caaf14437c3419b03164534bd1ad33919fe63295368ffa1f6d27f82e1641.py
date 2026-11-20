"""Preserve remaining historical shell commands; classify generated non-source files."""
from pathlib import Path
import ast,os,json,re,hashlib,datetime
ROOT=Path(__file__).resolve().parent
from circos_review_common import constrain_process
constrain_process()
PROJECT=Path(r'D:\Illumina\David_Eisele\terreus_genotyping')
OUT=ROOT/'repository_publication/aspergillus_terreus_genomics_pipeline'
versions=json.loads((OUT/'manifests/script_versions.json').read_text());syntax=json.loads((OUT/'manifests/historical_python_syntax.json').read_text())
objects={r['stored_sha256']:r['stored_path'] for r in versions};EXT={'.py','.sh','.bash','.R','.r','.ps1','.nf','.smk'}
parsed=ast.parse((ROOT/'prepare_git_publication.py').read_text())
exec(compile(ast.Module(body=[n for n in parsed.body if isinstance(n,ast.FunctionDef) and n.name in {'digest','category','archive'}],type_ignores=[]),'prepare_git_publication.py','exec'))
targets=[PROJECT/'inbox/ATSec/persvade.bash',PROJECT/'inbox/ksnp/unannotated_genomes/repeat_funannotat.bash']
targets += list(PROJECT.glob('results*/0_1_spades/*/*/run_spades.sh'))+list(PROJECT.glob('inbox/fastq/terreus/spades/*/run_spades.sh'))
assert len(targets)==71
known={(r['origin'],r['original_relative_path']) for r in versions}
for p in targets:
 if ('project',p.relative_to(PROJECT).as_posix()) not in known:archive(p,PROJECT,'project')
(OUT/'manifests/script_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
audit={'source_like_files_discovered_in_original_project':436,'genuine_extension_matched_script_files_preserved':83,'additional_extensionless_commands_not_in_that_count':'included in script manifest','excluded_generated_data_objects':{'df_benchmark_all.py':88,'df_cross_benchmark.py':88,'df_benchmark_all_with_integratedInfo.py':88,'output_rename-json.py':65,'reference_genome.fasta.chr_to_len.py':14,'rearranged_genome.fasta.windows81602bp.bed_GCcontent.py':2,'rearranged_genome.fasta.tmp.fasta.chr_to_len.py':2,'rearranged_genome.fasta.chr_to_len.py':2},'excluded_external_tool_copies':{'vcf2phylip.py':4},'evidence':'Generated dataframe .py sample begins with Python pickle protocol bytes (80 04), not Python source. PerSVade gene/chromosome/GC objects and annotation rename JSON caches are generated data, excluded under the no-data policy. External vcf2phylip copies are tool dependencies, not locally implemented code.','search_exclusions':['vendor','busco_downloads','.git','.snakemake','envs'],'scope':'All discovered genuine .py/.sh/.bash/.R/.r/.ps1/.nf/.smk source files in the audited project are accounted for; large generated result directories and bundled tool caches are not published.'}
(OUT/'manifests/original_project_source_audit.json').write_text(json.dumps(audit,indent=2))
s=json.loads((OUT/'manifests/publication_summary.json').read_text());s.update(script_source_records=len(versions),distinct_archived_source_objects=len(objects));(OUT/'manifests/publication_summary.json').write_text(json.dumps(s,indent=2))
print('source records',len(versions),'distinct objects',len(objects))
