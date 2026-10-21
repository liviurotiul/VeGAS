#!/usr/bin/env python3
"""Generate HTML phylogeny pages for every virus/reference group.

For each virus two pages are written: one with the reference sequence included
and one with only the samples.
"""
import argparse
import html
import os

# (subdirectory with tree.nwk, output page name, label)
VARIANTS = [
    ('with_reference', 'phylogeny_with_reference.html', 'probe + referință'),
    ('samples_only', 'phylogeny_samples_only.html', 'doar probe'),
]


def read_manifest(manifest_path):
    """Read the phylogeny manifest written by ``phylo.py``."""
    groups = []
    with open(manifest_path, 'rt') as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 4:
                continue
            reference, directory, sample_count, has_reference = parts
            groups.append((reference, directory, int(sample_count), has_reference == '1'))
    return groups


def inject_tree(template, newick):
    marker = '<script id="newick-data" type="text/plain"></script>'
    return template.replace(marker, f'<script id="newick-data" type="text/plain">{newick}</script>')


def index_html(groups, link_prefix=''):
    rows = []
    for reference, directory, sample_count, has_reference in groups:
        label = html.escape(reference)
        plural = 'probă' if sample_count == 1 else 'probe'
        links = []
        for subdir, page, text in VARIANTS:
            if subdir == 'with_reference' and not has_reference:
                continue
            href = html.escape(f'{link_prefix}{directory}/{page}')
            links.append(f'<a href="{href}">{text}</a>')
        rows.append(f'<li><strong>{label}</strong> ({sample_count} {plural}): ' + ' | '.join(links) + '</li>')
    content = '\n'.join(rows) if rows else '<li>Nu există probe cu o referință virală selectată.</li>'
    return f'''<!doctype html>
<html lang="ro"><head><meta charset="utf-8"><title>Arbori filogenetici</title>
<style>body{{font-family:system-ui,-apple-system,sans-serif;max-width:900px;margin:3rem auto;padding:0 1rem}} li{{margin:.7rem 0}}</style>
</head><body><h1>Arbori filogenetici</h1>
<p>Fiecare arbore conține exclusiv probele care au avut aceeași referință virală selectată la asamblare. </p>
<ul>{content}</ul></body></html>'''

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', default='data/output', help='Project output root')
    p.add_argument('--phylogeny-dir', default='phylogeny', help='Directory containing per-virus trees')
    p.add_argument('--template', default='src/html_templates/phylogeny.html', help='HTML template path')
    p.add_argument('--out', default='data/output/phylogeny.html', help='Output path for the phylogeny index page')
    args = p.parse_args()

    phylogeny_dir = os.path.join(args.root, args.phylogeny_dir)
    manifest_path = os.path.join(phylogeny_dir, 'trees.tsv')
    if not os.path.exists(manifest_path):
        print('Tree manifest not found:', manifest_path)
        return 1

    tmpl_path = args.template
    if not os.path.exists(tmpl_path):
        print('Template not found:', tmpl_path)
        return 2

    with open(tmpl_path, 'rt') as fh:
        tmpl = fh.read()

    groups = read_manifest(manifest_path)
    for reference, directory, _, _ in groups:
        for subdir, page, _ in VARIANTS:
            tree_path = os.path.join(phylogeny_dir, directory, subdir, 'tree.nwk')
            if not os.path.exists(tree_path):
                print('Tree file not found for', reference, f'({subdir}):', tree_path)
                continue
            with open(tree_path, 'rt') as fh:
                newick = fh.read().strip()
            output_path = os.path.join(phylogeny_dir, directory, page)
            with open(output_path, 'wt') as fh:
                fh.write(inject_tree(tmpl, newick))
            print('Wrote', output_path)

    index_path = os.path.join(phylogeny_dir, 'index.html')
    with open(index_path, 'wt') as fh:
        fh.write(index_html(groups))
    # The landing page lives outside phylogeny/, so its links need the folder prefix.
    out_dir = os.path.dirname(os.path.abspath(args.out))
    prefix = os.path.relpath(os.path.abspath(phylogeny_dir), out_dir).replace(os.sep, '/') + '/'
    with open(args.out, 'wt') as fh:
        fh.write(index_html(groups, prefix))
    print('Wrote', index_path)
    print('Wrote', args.out)

    return 0

if __name__ == '__main__':
    raise SystemExit(main())
