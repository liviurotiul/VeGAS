#!/usr/bin/env python3
"""Generate a simple static HTML page for the phylogeny tree.nwk

This script reads `phylogeny/tree.nwk` under the given root and writes
`data/output/phylogeny.html` (or custom out path) by injecting the Newick
string into the HTML template `src/html_templates/phylogeny.html`.

Usage:
  python scripts/generate_phylo_html.py --root data/output --out data/output/phylogeny.html
"""
import argparse
import os
import io

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', default='data/output', help='Project output root')
    p.add_argument('--tree', default='phylogeny/tree.nwk', help='Path to Newick tree inside root')
    p.add_argument('--template', default='src/html_templates/phylogeny.html', help='HTML template path')
    p.add_argument('--out', default='data/output/phylogeny.html', help='Output HTML path (base for generated files)')
    args = p.parse_args()

    tree_path = os.path.join(args.root, args.tree)
    if not os.path.exists(tree_path):
        print('Tree file not found:', tree_path)
        return 1

    tmpl_path = args.template
    if not os.path.exists(tmpl_path):
        print('Template not found:', tmpl_path)
        return 2

    with open(tree_path, 'rt') as fh:
        newick = fh.read().strip()

    # read template
    with open(tmpl_path, 'rt') as fh:
        tmpl = fh.read()

    # gather reference names (first token of FASTA header) from data/reference
    # prefer references located under the output root, but also check common repo locations
    ref_dir = os.path.join(args.root, 'reference')
    if not os.path.isdir(ref_dir):
        # try sibling 'reference' next to output root (e.g. data/reference)
        parent = os.path.dirname(args.root)
        if parent:
            alt = os.path.join(parent, 'reference')
            if os.path.isdir(alt):
                ref_dir = alt
    if not os.path.isdir(ref_dir):
        # final fallback: top-level data/reference
        alt2 = os.path.join('data', 'reference')
        if os.path.isdir(alt2):
            ref_dir = alt2
    ref_names = set()
    if os.path.isdir(ref_dir):
        for fn in os.listdir(ref_dir):
            if not fn.lower().endswith('.fasta') and not fn.lower().endswith('.fa'):
                continue
            fpath = os.path.join(ref_dir, fn)
            try:
                with open(fpath, 'rt') as fh:
                    for line in fh:
                        if line.startswith('>'):
                            # extract name up to whitespace
                            header = line[1:].strip().split()[0]
                            if header:
                                ref_names.add(header)
                            break
            except Exception:
                continue

    # simple Newick parser/serializer to allow filtering leaves by name
    class Node:
        def __init__(self):
            self.name = None
            self.length = None
            self.children = []

    def parse_newick(s):
        s = s.strip()
        i = 0
        L = len(s)
        def eat_ws():
            nonlocal i
            while i < L and s[i].isspace(): i += 1

        def parse_sub():
            nonlocal i
            eat_ws()
            node = Node()
            if i < L and s[i] == '(':
                i += 1
                while True:
                    child = parse_sub()
                    node.children.append(child)
                    eat_ws()
                    if i < L and s[i] == ',':
                        i += 1
                        continue
                    if i < L and s[i] == ')':
                        i += 1
                        break
                    break
            eat_ws()
            # parse name
            name = ''
            while i < L and s[i] not in ':,();\n\r\t':
                name += s[i]; i += 1
            if name:
                node.name = name
            eat_ws()
            if i < L and s[i] == ':':
                i += 1
                num = ''
                while i < L and s[i] in '0123456789.eE+-':
                    num += s[i]; i += 1
                try:
                    node.length = float(num) if num else None
                except Exception:
                    node.length = None
            return node

        root = parse_sub()
        return root

    def serialize(node):
        s = ''
        if node.children:
            s += '(' + ','.join(serialize(c) for c in node.children) + ')'
        if node.name:
            s += node.name
        if node.length is not None:
            s += ':' + ('%g' % node.length)
        return s

    def filter_node(node, refs):
        # returns node or None if filtered
        if not node.children:
            if node.name and node.name in refs:
                return None
            return node
        new_children = []
        for c in node.children:
            kept = filter_node(c, refs)
            if kept is not None:
                new_children.append(kept)
        node.children = new_children
        if not node.children:
            return None
        return node

    # prepare two Newick strings: original (with refs) and filtered (no refs)
    newick_clean = newick.rstrip(';')
    newick_with = newick_clean + ';'
    newick_no_refs = ''
    if ref_names:
        try:
            parsed = parse_newick(newick_clean)
            filtered = filter_node(parsed, ref_names)
            if filtered is None:
                newick_no_refs = ''
            else:
                newick_no_refs = serialize(filtered) + ';'
        except Exception:
            newick_no_refs = ''

    # helper to write one output file given a newick string
    def write_out(outpath, newick_str):
        if not newick_str:
            # write a page with a small message
            injected = tmpl.replace('<script id="newick-data" type="text/plain"></script>',
                                    '<script id="newick-data" type="text/plain"></script>')
            # add a small note in the container via JS: if no data, the existing script shows message
        else:
            injected = tmpl.replace('<script id="newick-data" type="text/plain"></script>',
                                    f'<script id="newick-data" type="text/plain">{newick_str}</script>')
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        with open(outpath, 'wt') as fh:
            fh.write(injected)
        print('Wrote', outpath)

    base_out = args.out
    base_dir = os.path.dirname(base_out) or args.root
    # produce two files
    out_with = os.path.join(base_dir, 'phylogeny_with_refs.html')
    out_no = os.path.join(base_dir, 'phylogeny_no_refs.html')
    write_out(out_with, newick_with)
    write_out(out_no, newick_no_refs)

    return 0

if __name__ == '__main__':
    raise SystemExit(main())
