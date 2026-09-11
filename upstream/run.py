#!/usr/bin/env python3
"""Isolated processed-data reproduction. Never executes inside the source checkout."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

HERE = Path(__file__).resolve().parent
BASE = '03_analysis_results/10_full_stable_id_rerun_20260823'
CROSS = '07_scripts/06_cross_tissue_integration/'
GO = BASE + '/03_mouse_go_concrete_terms_20260824/tables'
SEVEN = BASE + '/04_seven_pathway_response_26_tissues_20260824'
CAT = BASE + '/04_seven_pathway_modified_contract_updated_20260824'

def rows(name):
    with open(HERE / name, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def sha(p):
    with open(p, 'rb') as f:
        h = hashlib.sha256()
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
        return h.hexdigest()

def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)

def stages():
    return [
      ('01_go_mapping', 'Rscript', ['07_scripts/04_dna_damage_repair_focus/prepare_mouse_go_current_impl.R'], {}),
      ('02_go_concrete', 'Rscript', ['07_scripts/04_dna_damage_repair_focus/prepare_mouse_go_concrete_terms_current.R'], {}),
      ('03_stable_registry', 'Rscript', ['07_scripts/01_main_data_selection/build_stable_id_stage01_current.R'], {}),
      ('04_gsea', 'Rscript', ['07_scripts/02_main_global_analysis/run_stable_id_gsea_current.R'], {'RUN_GLOBAL_MSIGDB':'true','FGSEA_NPROC':'8'}),
      ('05_go_matrices', 'Rscript', ['07_scripts/02_main_global_analysis/summarize_explicit_GO_pathway_correlations_current.R'], {}),
      ('06_hallmark_selection', 'Rscript', ['07_scripts/02_main_global_analysis/plot_global_gsea_current.R'], {}),
      ('07_seven_pathway', 'Rscript', [CROSS+'run_seven_pathway_linear_matrices_current.R'], {}),
      ('08_seven_catalog', 'Rscript', [CROSS+'run_seven_pathway_linear_matrices_current.R','--outdir',CAT,'--rules','07_scripts/config/mouse_seven_pathway_go_seed_rules.csv','--dna-repair-relations','is_a,part_of'], {}),
      ('09_sample_counts', 'Rscript', [CROSS+'plot_all_tissue_seven_pathway_member_counts_per_sample_20260829.R'], {}),
      ('10_dsb_meta', 'Rscript', [CROSS+'plot_essential_components_counts_logfc_matrix.R'], {}),
      ('11_all_component_meta', 'Rscript', [CROSS+'plot_essential_components_log2fc_pvalue_heatmaps_current.R'], {}),
      ('12_ssb_meta', 'Rscript', [CROSS+'plot_ssb_meta_log2fc_without_neil2_20260829.R'], {}),
      ('13_qsmooth', 'Rscript', [CROSS+'plot_essential_components_yarn_qsmooth_expression_heatmaps_current.R'], {}),
      ('14_remove_neil2', 'Rscript', [CROSS+'rebuild_essential_components_without_neil2_20260828.R'], {}),
      ('15_ddr_membership', 'python', [CROSS+'build_mouse_go_ddr.py'], {
          'MOUSE_GO_DDR_OUTPUT_DIR':BASE+'/05_direct_ddr_final_fdr',
          'MOUSE_GO_MAPPING_FILE':'08_GO_annotation/mgi_mod_gaf_five_relation_concrete_terms_20260824/final_manual_curated/03_mouse_GO_annotations_with_gene_ids.tsv.gz',
          'MOUSE_GO_ANALYSIS_VARIANT':'mgi_mod_gaf_five_relation_concrete_terms_20260824',
          'MOUSE_GO_STABLE_LONG_STATS':BASE+'/01_stable_id_registry/03_full_gene_statistics_ensembl_long.csv.gz',
          'MOUSE_GO_EXACT_DATASET_RESULTS':BASE+'/01_stable_id_registry/01_dataset_manifest_validated.csv'}),
      ('16_common_direction', 'python', [CROSS+'build_ddr_common_direction_four_groups_log2fc_current.py'], {}),
      ('17_venn', 'python', [CROSS+'build_ddr_go0006974_dataset_level_venn_20260826.py'], {'DDR_VENN_OUTPUT_DIR':'03_analysis_results/19_ddr_go0006974_dataset_level_venn_20260826'}),
      ('17b_upset_tables', 'Rscript', [CROSS+'plot_teacher_selected_go_upsets_20260829.R'], {}),
      ('18_mouse_metadata', 'Rscript', ['publication_release/R/data/01_refresh_osdr_mouse_metadata.R','--scope={work}/03_analysis_results/24_teacher_selected_figure_release_20260829/seven_pathway_counts/tables/02_seven_pathway_member_gene_counts_per_source_sample_long.csv'], {}),
      ('19_metadata_audit', 'Rscript', ['publication_release/R/data/02_build_mouse_sample_metadata_audit.R'], {}),
      ('20_sample_spearman', 'python', [CROSS+'build_log2fc_per_tissue_pathway_spearman_20260902.py'], {}),
    ]

def init(work):
    if work.exists():
        raise SystemExit('init requires a new, nonexistent directory; select a new --work path')
    shutil.copytree(HERE/'sources', work)
    for p in (HERE/'design_inputs').rglob('*'):
        if p.is_file():
            copy(p, work/p.relative_to(HERE/'design_inputs'))
    p=work/'03_analysis_results/16_essential_components_yarn_qsmooth_expression_heatmaps_20260826/input/Essential_components_1_teacher_source.xlsx'
    copy(p,work/'design/Essential_components_1_teacher_source.xlsx')
    (work/'.upstream-reproduction').write_text('isolated processed-data reconstruction\n')
    (work/'logs').mkdir()
    print('Initialized', work)

def check_work(work):
    if not (work/'.upstream-reproduction').is_file():
        raise SystemExit('Not an initialized reproduction workspace. Run init in a NEW directory.')
    if any(p.is_symlink() for p in work.rglob('*')
           if p.relative_to(work).parts[0] not in {'.venv', '.runtime', 'R'}):
        raise SystemExit('Reproduction input/output trees must not contain symlinks')

def acquire(work, local, download):
    archive = HERE.parent / "data/annotation_snapshot_20260822.tar.gz"
    if archive.is_file():
        expected={r['workspace_path']:r for r in rows('annotation_snapshots.csv')}
        with tarfile.open(archive, 'r:gz') as bundle:
            members={m.name.removeprefix('./'):m for m in bundle.getmembers() if not m.isdir()}
            if set(members)!=set(expected):
                raise SystemExit('Annotation archive member inventory differs from manifest')
            for rel, row in expected.items():
                member=members[rel]
                if not member.isfile() or Path(rel).is_absolute() or '..' in Path(rel).parts:
                    raise SystemExit('Unsafe annotation archive member: '+rel)
                dst=work/rel
                if dst.exists():
                    if sha(dst)!=row['sha256']:raise SystemExit('Existing snapshot changed: '+rel)
                    continue
                dst.parent.mkdir(parents=True,exist_ok=True)
                partial=dst.with_name(dst.name+'.partial')
                with bundle.extractfile(member) as src,open(partial,'wb') as out:
                    shutil.copyfileobj(src,out)
                if sha(partial)!=row['sha256']:
                    partial.unlink()
                    raise SystemExit('Annotation archive checksum mismatch: '+rel)
                partial.replace(dst)
    errors=[]
    for row in rows('nasa_processed_files.csv')+rows('annotation_snapshots.csv'):
        dst=work/row['workspace_path']
        if dst.is_file() and sha(dst)==row['sha256']:
            continue
        if dst.exists():
            errors.append('Existing file checksum differs: '+str(dst));continue
        source=local/row['workspace_path'] if local else None
        if source and source.is_file():
            if sha(source)!=row['sha256']:
                errors.append('Source checksum differs: '+str(source));continue
            copy(source,dst);continue
        if download and row.get('files_api'):
            try:
                req=urllib.request.Request(row['files_api'],headers={'User-Agent':'publication-reproduction/1.0'})
                with urllib.request.urlopen(req,timeout=90) as res:
                    data=json.load(res)
                url=data[row['accession']]['files'][row['filename']]['URL']
                dst.parent.mkdir(parents=True,exist_ok=True)
                partial=dst.with_name(dst.name+'.partial')
                with urllib.request.urlopen(url,timeout=180) as res,open(partial,'wb') as f:
                    shutil.copyfileobj(res,f)
                if sha(partial)!=row['sha256']:
                    partial.unlink()
                    raise ValueError('NASA file changed: downloaded SHA-256 differs from approved input')
                partial.replace(dst)
                with open(work/'logs/acquisition.jsonl','a') as f:
                    f.write(json.dumps({'path':row['workspace_path'],'url':url,'sha256':row['sha256']})+'\n')
            except Exception as e:
                errors.append(row['workspace_path']+': '+str(e))
        else:
            errors.append('Exact archival file required: '+row['workspace_path'])
    mapping='ensembl116_full_mgi_entrez_mapping.tsv.gz'
    src=work/'08_GO_annotation/mgi_mod_gaf_five_relation_20260822/source'/mapping
    if src.exists():copy(src,work/'08_GO_annotation/mgi_mod_gaf_20260804'/mapping)
    for row in rows('msigdb_snapshot.csv'):
        dst=work/row['workspace_path']
        options=[HERE.parent/'data'/Path(row['workspace_path']).name]
        if local:options.insert(0,local/row['local_source_path'])
        if not dst.exists():
            source=next((p for p in options if p.is_file()),None)
            if source:
                if sha(source)!=row['sha256']:raise SystemExit('MSigDB snapshot checksum mismatch')
                copy(source,dst)
        if dst.exists() and sha(dst)!=row['sha256']:raise SystemExit('MSigDB snapshot checksum mismatch')
    if errors:raise SystemExit('\n'.join(errors))
    print('All NASA processed files and historical annotation snapshots verified.')

def bridge_for_spearman(work):
    for n in ['04_tissue_statistics_concrete_terms_and_context.csv','13_tissue_pathway_Spearman_rho_matrix.csv']:
        copy(work/GO/n,work/'publication_release/data/publication_input/go'/n)

def run(work, first, last, rscript, python, library):
    seq=stages();names=[s[0] for s in seq]
    start=names.index(first) if first else 0;end=names.index(last)+1 if last else len(seq)
    if end<=start:raise SystemExit('--until precedes --from-stage')
    for name,program,args,overrides in seq[start:end]:
        if name=='20_sample_spearman':bridge_for_spearman(work)
        # Remove inherited historical output/root overrides before invoking originals.
        env={k:v for k,v in os.environ.items() if not k.startswith(('MOUSE_GO_','STABLE_','DDR_','GO_TERM_','THYMUS_','PROLIFERATION_','SEVEN_PATHWAY_','TEACHER_RELEASE_','SSB_META_','PUBLICATION_OUTPUT_'))}
        env.update(GRADUATE_DESIGN_ROOT=str(work),PLOT_GO_PATHWAY_LINEAR_MATRIX_SKILL=str(work/'vendor/plot-go-pathway-linear-matrix'),PYTHON=python,PYTHONPATH=str(work/'07_scripts'),PYTHONHASHSEED='25',MPLCONFIGDIR=str(work/'.matplotlib'),TMPDIR=tempfile.gettempdir())
        (work/'tmp').mkdir(exist_ok=True)
        if library:env['R_LIBS_USER']=library
        env.update(overrides)
        membership=work/'reference/MSigDB_2026.1.Mm_mouse_Ensembl_membership.tsv.gz'
        if membership.exists():env['MSIGDB_MEMBERSHIP_FILE']=str(membership)
        command=[rscript if program=='Rscript' else python]+[a.replace('{work}',str(work)) for a in args]
        print(name,flush=True)
        with open(work/'logs'/f'{name}.log','w') as log:
            result=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            raise SystemExit(f'{name} failed ({result.returncode}); see {work}/logs/{name}.log')
        (work/'logs'/f'{name}.done').write_text(json.dumps(command)+'\n')
    print('Selected original computation stages completed. Run export; then compare candidate inputs.')

def export(work):
    out=work/'candidate_publication';out.mkdir(exist_ok=True)
    report=[]
    for row in rows('approved_artifact_map.csv'):
        matches=[work/p for p in row['source_candidates'].split(';') if p and (work/p).is_file()]
        src=next((p for p in matches if sha(p)==row['sha256']),matches[0] if matches else None)
        if src:
            copy(src,out/row['publication_path'])
        report.append({'publication_path':row['publication_path'],'expected_sha256':row['sha256'],'candidate_sha256':sha(src) if src else '', 'source':str(src.relative_to(work)) if src else '', 'status':'byte_identical' if src and sha(src)==row['sha256'] else 'needs_semantic_comparison' if src else 'missing'})
    metadata=work/'publication_release/data/publication_input/mouse_metadata'
    if metadata.exists():shutil.copytree(metadata,out/'data/publication_input/mouse_metadata',dirs_exist_ok=True)
    audit=work/'publication_release/results/tables/Fig_6_mouse_sample_metadata_complete_audit.csv'
    if audit.exists():copy(audit,out/'results/tables'/audit.name)
    with open(work/'logs/export_comparison.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(report[0]));w.writeheader();w.writerows(report)
    print('Candidate only:',out)
    print({s:sum(r['status']==s for r in report) for s in {r['status'] for r in report}})
    for rule in rows('approved_label_overrides.csv'):
        target=out/rule['publication_path']
        if not target.exists():continue
        with target.open(newline='') as f:
            content=list(csv.DictReader(f)); fields=list(content[0])
        for record in content:
            if record[rule['key_column']]==rule['key_value'] and record[rule['column']]==rule['from_value']:
                record[rule['column']]=rule['to_value']
        with target.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(content)
    if any(r['status']=='missing' for r in report):
        raise SystemExit('Missing computed publication inputs; see logs/export_comparison.csv')

def render(work, rscript, python, library):
    """Render solely from exported computed inputs, plus presentation code."""
    export(work)
    candidate=work/'candidate_publication'
    for directory in ['R','python','workflow']:
        shutil.copytree(HERE.parent/directory,candidate/directory,dirs_exist_ok=True)
    shutil.copytree(HERE.parent/'data/external_statistics',candidate/'data/external_statistics',dirs_exist_ok=True)
    copy(HERE.parent/'provenance/final_figure_sources.csv',candidate/'provenance/final_figure_sources.csv')
    copy(HERE.parent/'README.md',candidate/'README.md')
    env=os.environ.copy();env['PYTHON']=python
    if library:env['R_LIBS_USER']=library
    if Path(rscript).is_absolute():env['PATH']=str(Path(rscript).parent)+os.pathsep+env.get('PATH','')
    subprocess.run([rscript,str(candidate/'workflow/run_final_release.R'),str(work/'reproduced_figures')],env=env,check=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['plan','init','acquire','run','export','render'])
    p.add_argument('--work',type=Path,default=HERE/'work')
    p.add_argument('--local-source',type=Path)
    p.add_argument('--download',action='store_true')
    p.add_argument('--from-stage');p.add_argument('--until')
    p.add_argument('--rscript',default='Rscript');p.add_argument('--python',default=sys.executable)
    p.add_argument('--r-library',help='Existing R package library; no packages installed automatically')
    a=p.parse_args();work=a.work.resolve()
    if a.action=='plan':
        for name,prog,args,env in stages():print(name,prog,*args,json.dumps(env))
    elif a.action=='init':init(work)
    else:
        check_work(work)
        if a.action=='acquire':acquire(work,a.local_source.resolve() if a.local_source else None,a.download)
        elif a.action=='run':run(work,a.from_stage,a.until,a.rscript,a.python,a.r_library)
        elif a.action=='export':export(work)
        elif a.action=='render':render(work,a.rscript,a.python,a.r_library)

if __name__=='__main__':main()
