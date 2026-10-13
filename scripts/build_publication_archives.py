#!/usr/bin/env python3
"""Build one autonomous publication archive per article; never publish it.

Uses standard-library Python. --verify runs the stated bounded CPU routes in
freshly extracted archives, without a GPU or a new large character scan.
The complete atlas release (4.85 GB) is an external data component, identified
by its original manifest. --legacy-archive bundles the original v1.0 deposit
unchanged in the certificates package when that ZIP is locally available.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ARTICLES = [
    ('certificates', 'RABUNG_BOUNDS_EN', 'v2.0', 'rabung_certificates_v2_publication'),
    ('gpu', 'GPU_METHOD_EN', 'v1.0', 'rabung_gpu_method_v1_publication'),
    ('atlas', 'PAPER_EN', 'v1.0', 'rabung_atlas_v1_publication'),
    ('cubic', 'CUBIC_QUADRUPLETS_EN', 'v1.0', 'cubic_quadruplets_v1_publication'),
]
COMMON = ['RIGHTS.md', 'docs/THIRD_PARTY.md', 'scripts/build_publication_archives.py']
COMMANDS = {
    'certificates': [['scripts/verify_bounds_update.py']],
    'gpu': [['scripts/verify_gpu_method.py']],
    'atlas': [['scripts/verify_atlas.py'], ['scripts/verify_ternary.py'],
              ['reproducibility/proofs/verify_synthesis.py'],
              ['reproducibility/proofs/verify_report4.py'],
              ['reproducibility/anteriorities/theorie_densites_residuels_20261004/preuves/verify_results.py'],
              ['scripts/verify_editorial.py', '--scope', 'atlas', '--catalogue']],
    'cubic': [['scripts/verify_editorial.py', '--scope', 'cubic'],
              ['scripts/cubic_quadruplets.py', '1000033']],
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def all_files(directory):
    result=[]
    for p in sorted((ROOT/directory).rglob('*')):
        if p.is_file() and not p.is_symlink() and '__pycache__' not in p.parts:
            result.append(p.relative_to(ROOT).as_posix())
    if not result:
        raise FileNotFoundError(directory)
    return result

def selection(kind, stem):
    paths=COMMON+[f'paper/{stem}.tex',f'paper/{stem}.pdf']
    if kind=='certificates':
        paths+=all_files('reproducibility/bounds_update')
        paths+=['scripts/verify_bounds_update.py','data/derived/ATLAS_CATALOGUE_K3_K64.csv']
    elif kind=='gpu':
        paths+=all_files('reproducibility/gpu_method')+['scripts/verify_gpu_method.py']
    elif kind=='cubic':
        paths+=all_files('reproducibility/cubic_quadruplets')
        paths+=['scripts/verify_editorial.py','scripts/verify_hypergeometric.py','scripts/cubic_quadruplets.py']
    elif kind=='atlas':
        for directory in ['analysis','src','data/derived','data/sample','data/ternary-million',
                          'reproducibility/anteriorities','reproducibility/evidence/new',
                          'reproducibility/proofs','reproducibility/ternary_diagnostics']:
            paths+=all_files(directory)
        paths+=['data/atlas/manifest.json','docs/DATA_FORMAT.md','docs/METHODS.md',
                'docs/REPRODUCTION.md','docs/CODE_EXTRACTION.json']
        paths+=[f'scripts/{name}.py' for name in [
            'atlas_codec','atlas_io','build_ternary_worker','build_workers','download_data',
            'query','query_ternary','recompute_small','reproduce_analysis','ternary_capture',
            'ternary_graph','ternary_io','verify_atlas','verify_ternary','verify_editorial']]
        paths+=all_files('tests')
    if kind in ('atlas','cubic'):
        paths+=['docs/EDITORIAL_PROVENANCE.json','docs/EDITORIAL_SHA256.json']
    return sorted(set(paths))

def inline_abstract(source):
    body=re.search(r'\\begin\{abstract\}([\s\S]*?)\\end\{abstract\}',source).group(1)
    body=' '.join(body.split())
    body=re.sub(r'\$([^$]+)\$',lambda m:r'\('+m.group(1)+r'\)',body)
    return body.replace('--','-')+'\n'

def readme(kind, stem, version, source_commit, legacy):
    limits={
        'certificates':'The default checks inspect retained complete CPU records and local powers; they do not repeat the billion-position scans. The original thirteen-claim gate and closure remain unchanged.',
        'gpu':'The default command reproduces the two timing tables, retained profiles and finite source bridges. It does not compile or execute CUDA and does not reproduce the timing experiment. The selected packet identifies additional upstream sources explicitly.',
        'atlas':'The complete atlas data comprise 4,852,996,533 bytes in 72 separate release assets. They are not embedded in this article archive. data/atlas/manifest.json preserves every asset size and SHA256, and scripts/download_data.py is the download route. The full dataset is publicly available from https://github.com/Vulkin-prog/rabung-atlas32/releases/tag/atlas32-v1.0.0. The one million ternary profiles, their retained positions and the distinct p < 2^34 collector are included here. The two sample atlas blocks are non-contiguous examples, not a smaller census.',
        'cubic':'The package needs no atlas data or GPU. Its programs verify finite identities, trace certificates and the low-degree moment bookkeeping. They do not computationally prove the cited automorphy or prime-number theorems. Four moments are not a full limiting distribution or a law of logarithmic-length maxima.',
    }
    lines=[f'# {stem}: publication archive {version}','',
           f'Author: Brice Pouly. Prepared edition: 12 October 2026.',
           f'Source repository snapshot: `{source_commit}`. This is an archive build, not a publication receipt.','',
           f'The PDF and standalone LaTeX source are in `paper/`. `ABSTRACT_ZENODO.txt` contains the abstract with inline LaTeX delimiters. `MANIFEST_SHA256.json` binds all included files.','',
           '## Reproduce the article','',
           'Run from this archive root. Use ordinary Python 3.10+; do not use `python -O` for inherited assertion-based checks.','',
           '```bash',*['python3 -B '+' '.join(args) for args in COMMANDS[kind]],'```','',
           ('The catalogue replay needs a C++17 compiler. Bounded CPU rebuild commands are in `docs/REPRODUCTION.md`.' if kind=='atlas' else None), '',
           'To rebuild the PDF with TeX Live (including lmodern, tikz and pgfplots where used):','',
           '```bash','mkdir -p build/paper',
           f'pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/paper paper/{stem}.tex',
           f'pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/paper paper/{stem}.tex','```','',
           '## Scope and evidence','',limits[kind],'',
           'Paths in the article refer to this layout unless an original historical package or upstream snapshot is explicitly named. Internal review reports and the other articles in the series are excluded.','',
           'For atlas/cubic packages, `docs/EDITORIAL_SHA256.json` is the scoped subset of the repository manifest for files included here. The full byte manifest is `MANIFEST_SHA256.json`.','',
           '## Rights and citations','',
           'Existing component notices are retained. RIGHTS.md records the current rights of newly assembled material; this packaging applies no new licence. The certificates paper has published editions v1.0.0 (DOI 10.5281/zenodo.22980965) and v2.0 (DOI 10.5281/zenodo.23337824). The GPU method paper is published as v1.0 (DOI 10.5281/zenodo.23338707). Series references identify the edition used; further editions are added only after publication.','']
    if kind=='certificates':
        lines+=['The nine legacy inputs in `reproducibility/bounds_update/published_v1_0/` are the exact copied inputs needed for the closure check.']
        lines+=[('The complete original deposited ZIP is included unchanged as `historical-publication/vdw_rabung_publication_package_v1.zip`. Extract it separately to replay the old campaign and empirical appendices; its original paths and component licence texts are preserved.' if legacy else 'The complete original campaign and empirical appendix artifact remains separately available at the published v1.0 DOI; it is not reconstructed from the nine closure inputs.'),'']
    return '\n'.join(x for x in lines if x is not None)

def verify_integrity(root):
    records=json.loads((root/'MANIFEST_SHA256.json').read_text())
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    if actual != set(records)|{'MANIFEST_SHA256.json'}:
        raise ArithmeticError('Archive file inventory differs')
    for name,item in records.items():
        p=root/name
        if p.stat().st_size!=item['bytes'] or sha(p)!=item['sha256']:
            raise ArithmeticError('Archive integrity differs: '+name)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'build/publication')
    parser.add_argument('--pdf-dir',type=Path,default=ROOT/'paper')
    parser.add_argument('--legacy-archive',type=Path)
    parser.add_argument('--source-commit',default=os.environ.get('GITHUB_SHA','4700f6673d468c4767846b9e433d740f640d35cb (editorial baseline; see archive metadata)'))
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    for kind,stem,version,name in ARTICLES:
        with tempfile.TemporaryDirectory(prefix='rabung-publication-') as temp:
            folder=Path(temp)/name;folder.mkdir()
            paths=selection(kind,stem)
            for relative in paths:
                src=args.pdf_dir/f'{stem}.pdf' if relative==f'paper/{stem}.pdf' else ROOT/relative
                if not src.is_file() or src.is_symlink():raise FileNotFoundError(src)
                dst=folder/relative;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
            if kind in ('atlas','cubic'):
                mp=folder/'docs/EDITORIAL_SHA256.json'
                original=json.loads(mp.read_text())
                mp.write_text(json.dumps({p:sha(folder/p) for p in original if (folder/p).is_file()},indent=2,sort_keys=True)+'\n')
            if kind=='certificates' and args.legacy_archive:
                dst=folder/'historical-publication/vdw_rabung_publication_package_v1.zip';dst.parent.mkdir()
                shutil.copyfile(args.legacy_archive,dst)
            (folder/'ABSTRACT_ZENODO.txt').write_text(inline_abstract((folder/f'paper/{stem}.tex').read_text()))
            if kind=='certificates':
                (folder/'VERSION_NOTE_ZENODO.txt').write_text('Version 2 adds six primitive certificates, rechecked by complete CPU scans, yielding fourteen direct bounds and stronger three-colour and six-colour recurrence consequences. It revises the exposition and reproducibility package while preserving the original thirteen-claim evidence and closed empirical analyses.\n')
            (folder/'README.md').write_text(readme(kind,stem,version,args.source_commit,bool(args.legacy_archive)))
            meta={'article':kind,'version':version,'date':'2026-10-12','source_repository':'Vulkin-prog/rabung-atlas32',
                  'source_commit':args.source_commit,'source_sha256':sha(folder/f'paper/{stem}.tex'),
                  'pdf_sha256':sha(folder/f'paper/{stem}.pdf'),'publication_status':'prepared; not published by this build',
                  'standalone_latex':True,'requires_other_unpublished_articles':False,
                  'external_complete_atlas':kind=='atlas','legacy_deposit_included':kind=='certificates' and bool(args.legacy_archive)}
            (folder/'PUBLICATION_METADATA.json').write_text(json.dumps(meta,indent=2)+'\n')
            records={p.relative_to(folder).as_posix():{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(folder.rglob('*')) if p.is_file()}
            (folder/'MANIFEST_SHA256.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')
            target=out/f'{name}.zip'
            with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
                for p in sorted(folder.rglob('*')):
                    if p.is_file():
                        info=zipfile.ZipInfo(p.relative_to(Path(temp)).as_posix(),date_time=(2026,10,12,12,0,0))
                        info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
                        z.writestr(info,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=6)
            with tempfile.TemporaryDirectory(prefix='rabung-unpacked-check-') as unpack:
                with zipfile.ZipFile(target) as z:
                    if z.testzip() is not None:raise ArithmeticError('ZIP CRC failure')
                    z.extractall(unpack)
                extracted=Path(unpack)/name;verify_integrity(extracted)
                if args.verify:
                    for command in COMMANDS[kind]:
                        subprocess.run([sys.executable,'-B',*command],cwd=extracted,check=True,timeout=900)
            print(json.dumps({'article':kind,'archive':target.name,'files':len(records)+1,'bytes':target.stat().st_size,
                              'sha256':sha(target),'integrity':'PASS','bounded_commands_replayed':args.verify}),flush=True)

if __name__=='__main__':main()
