#!/usr/bin/env python3
"""Generate an IGV HTML report listing all assembly BAMs.

Scans for BAM files (defaults to data/output/assembly), builds a `rows` list with
relative URLs (suitable when serving `data/output` as the web root) and renders
an HTML file from the template `src/html_templates/igv.html` by replacing the
placeholder marker `__ROWS_JSON__` with the generated JSON.

Usage:
  python src/scripts/igv_report.py --root data/output --out data/output/igv_report.html
"""
import os
import glob
import json
import argparse
from pathlib import Path


def find_bams(root_dirs):
    paths = []
    for d in root_dirs:
        paths.extend(glob.glob(os.path.join(d, "*.bam")))
    # dedupe and sort
    filtered = []
    for path in sorted(set(paths)):
        if "biolm" in path.lower():
            continue
        filtered.append(path)
    return filtered


def rel_to_root(path, root):
    try:
        return os.path.relpath(path, root).replace('\\', '/')
    except Exception:
        return os.path.basename(path)


def guess_fai_for(fasta_path):
    if not fasta_path:
        return None
    p = Path(fasta_path)
    maybe = str(p) + ".fai"
    if os.path.exists(maybe):
        return maybe
    # sometimes .fai next to fasta but different extension; fallback None
    return None


def first_ref_from_fai(fai_path):
    if not fai_path or not os.path.exists(fai_path):
        return None
    with open(fai_path, 'rt') as fh:
        first = fh.readline().strip().split('\t')[0]
        return first


def build_rows(bam_files, web_root):
    rows = []
    for bam in bam_files:
        base = os.path.splitext(os.path.basename(bam))[0]
        if "biolm" in base.lower():
            continue
        bai_candidates = [bam + '.bai', os.path.splitext(bam)[0] + '.bai']
        bai = next((b for b in bai_candidates if os.path.exists(b)), None)

        # try to find a fasta next to bam (assembly/sample.fasta) or in ../reference
        fasta_candidates = [
            os.path.join(os.path.dirname(bam), base + '.fasta'),
            os.path.join(os.path.dirname(bam), base + '.fa'),
            os.path.join('data', 'reference', base + '.fasta'),
            os.path.join('data', 'reference', base + '.fa'),
        ]
        fasta = next((f for f in fasta_candidates if os.path.exists(f)), None)
        fai = guess_fai_for(fasta) if fasta else None

        first_ref = first_ref_from_fai(fai) if fai else None

        row = {
            'name': base,
            'bam': rel_to_root(bam, web_root),
            'bai': rel_to_root(bai, web_root) if bai else None,
            'fasta': rel_to_root(fasta, web_root) if fasta else None,
            'fai': rel_to_root(fai, web_root) if fai else None,
            'locus': first_ref or None
        }
        rows.append(row)
    return rows


def render_template(template_path, rows_json, out_path):
    with open(template_path, 'rt', encoding='utf-8') as fh:
        tpl = fh.read()
    if '__ROWS_JSON__' not in tpl:
        raise RuntimeError('Template must contain the placeholder __ROWS_JSON__')
    out_html = tpl.replace('__ROWS_JSON__', rows_json)
    # Ensure the output directory exists. If out_path has no dirname (it's
    # intended to be created in the current working directory), use '.' so
    # os.makedirs doesn't receive an empty string which raises on some OSes.
    out_dir = os.path.dirname(out_path) or '.'
    os.makedirs(out_dir, exist_ok=True)
    with open(out_path, 'wt', encoding='utf-8') as fh:
        fh.write(out_html)


def main():
    p = argparse.ArgumentParser(description='Generate IGV report HTML from BAMs')
    p.add_argument('--root', default='data/output', help='Web root (where report will be served from)')
    p.add_argument('--template', default='src/html_templates/igv.html', help='HTML template path')
    p.add_argument('--out', default='data/output/igv_report.html', help='Output HTML path')
    p.add_argument('--bamdirs', nargs='*', default=['data/output/assembly', 'data/output'], help='Directories to scan for BAM files')
    args = p.parse_args()

    # Resolve bamdirs: accept absolute paths, existing relative paths, or
    # interpret relative paths as relative to the provided --root.
    resolved_dirs = []
    for d in args.bamdirs:
        if os.path.isabs(d):
            resolved_dirs.append(d)
        elif os.path.exists(d):
            resolved_dirs.append(os.path.abspath(d))
        else:
            # interpret relative to root
            candidate = os.path.join(args.root, d)
            resolved_dirs.append(candidate)

    bam_files = find_bams(resolved_dirs)
    if not bam_files:
        print('No BAM files found in', args.bamdirs)
    rows = build_rows(bam_files, args.root)
    rows_json = json.dumps(rows, indent=2)
    render_template(args.template, rows_json, args.out)
    print(f'Wrote IGV report to {args.out} with {len(rows)} rows')


if __name__ == '__main__':
    main()
