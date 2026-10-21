#!/usr/bin/env python3
"""Build one phylogeny per virus/reference selected during assembly.

Each sample is assigned to the reference written in
``assembly/<sample>.best_ref.txt``. Only samples with the same selected
reference are aligned together, so sequences from different viruses are never
combined in one tree. The output directory contains one subdirectory per
reference and a ``trees.tsv`` manifest for the HTML report.

For every virus two trees are built:
  - ``with_reference/``: the samples plus the reference sequence (leaf ``REF_<name>``)
  - ``samples_only/``: only the samples
"""
import argparse
import glob
import os
import re
import subprocess


def read_first_record(fasta_path):
    """Return tuple (header, seq) of the first FASTA record in file, or None."""
    header = None
    parts = []
    with open(fasta_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                header = line[1:].split()[0]
                break
        if header is None:
            return None
        # read rest of first seq
        for line in fh:
            if line.startswith('>'):
                break
            parts.append(line.strip())
    seq = ''.join(parts)
    return header, seq


def write_multi_fasta(records, out_path):
    with open(out_path, 'wt') as fh:
        for name, seq in records:
            fh.write(f'>{name}\n')
            # wrap at 80
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')


def run_cmd(cmd, cwd=None):
    print('Running:', ' '.join(cmd))
    res = subprocess.run(cmd, cwd=cwd)
    if res.returncode != 0:
        raise RuntimeError(f'Command failed: {cmd}')


def reference_name(best_ref_path):
    """Return the reference basename recorded by the assembly step."""
    with open(best_ref_path, 'rt') as fh:
        reference = fh.read().strip()
    if not reference:
        raise ValueError('empty best-reference file')
    return os.path.basename(reference).removesuffix('.fasta')


def safe_directory_name(name):
    """Produce a portable directory name while preserving readable virus names."""
    return re.sub(r'[^A-Za-z0-9._-]+', '_', name).strip('._') or 'unknown_reference'


def read_all_records(fasta_path):
    """Return a list of (name, sequence) for every record in a FASTA file."""
    records = []
    name = None
    parts = []
    with open(fasta_path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(parts)))
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        records.append((name, ''.join(parts)))
    return records


def find_reference_sequence(reference, best_ref_value, reference_dir):
    """Return the first sequence of the reference FASTA, or None if not found."""
    candidates = [best_ref_value + '.fasta', os.path.join(reference_dir, reference + '.fasta')]
    for path in candidates:
        if os.path.exists(path):
            rec = read_first_record(path)
            if rec is not None:
                return rec[1]
    return None


def pairwise_newick(alignment_path):
    """Newick for exactly two aligned sequences (FastTree needs more leaves)."""
    (name_a, seq_a), (name_b, seq_b) = read_all_records(alignment_path)
    compared = [(x, y) for x, y in zip(seq_a.upper(), seq_b.upper()) if x != '-' and y != '-']
    diffs = sum(1 for x, y in compared if x != y)
    half = (diffs / len(compared) / 2) if compared else 0.0
    return f'({name_a}:{half:.6f},{name_b}:{half:.6f});\n'


def build_tree(records, out_dir, mafft_opts):
    """Write inputs, alignment and Newick tree for one group of sequences."""
    os.makedirs(out_dir, exist_ok=True)
    inputs_path = os.path.join(out_dir, 'inputs.fasta')
    alignment_path = os.path.join(out_dir, 'alignment.fa')
    tree_path = os.path.join(out_dir, 'tree.nwk')
    write_multi_fasta(records, inputs_path)

    if len(records) == 1:
        # A single sequence has no topology; keep it as a one-leaf tree.
        write_multi_fasta(records, alignment_path)
        with open(tree_path, 'wt') as fh:
            fh.write(f'{records[0][0]};\n')
        return tree_path

    with open(alignment_path, 'wt') as fh:
        result = subprocess.run(['mafft'] + mafft_opts.split() + [inputs_path], stdout=fh)
    if result.returncode != 0:
        raise RuntimeError('MAFFT failed with return code ' + str(result.returncode))

    if len(records) == 2:
        with open(tree_path, 'wt') as fh:
            fh.write(pairwise_newick(alignment_path))
        return tree_path

    with open(tree_path, 'wt') as fh:
        result = subprocess.run(['FastTree', '-nt', alignment_path], stdout=fh)
    if result.returncode != 0:
        with open(tree_path, 'wt') as fh:
            result = subprocess.run(['fasttree', '-nt', alignment_path], stdout=fh)
        if result.returncode != 0:
            raise RuntimeError('FastTree failed')
    return tree_path


def main():
    p = argparse.ArgumentParser(description='Create one phylogenetic tree per selected reference')
    p.add_argument('--assembly-dir', default='assembly', help='Directory with per-sample FASTA files')
    p.add_argument('--best-ref-dir', default='assembly', help='Directory with per-sample .best_ref.txt files')
    p.add_argument('--reference-dir', default='ref_indexes', help='Directory with reference <name>.fasta files')
    p.add_argument('--out-dir', default='phylogeny', help='Output directory for alignment and tree')
    p.add_argument('--mafft-opts', default='--auto', help='Options to pass to MAFFT')
    args = p.parse_args()

    fasta_files = sorted(glob.glob(os.path.join(args.assembly_dir, '*.fasta')))
    if not fasta_files:
        print('No assembly fasta files found in', args.assembly_dir)
        return

    records_by_reference = {}
    best_ref_values = {}
    for f in fasta_files:
        sample = os.path.splitext(os.path.basename(f))[0]
        rec = read_first_record(f)
        if rec is None:
            print('Warning: no FASTA records in', f)
            continue
        best_ref_path = os.path.join(args.best_ref_dir, f'{sample}.best_ref.txt')
        if not os.path.exists(best_ref_path):
            print('Warning: no selected reference for', sample, '- skipping phylogeny')
            continue
        try:
            reference = reference_name(best_ref_path)
        except ValueError as error:
            print(f'Warning: invalid selected reference for {sample}: {error}')
            continue
        with open(best_ref_path, 'rt') as fh:
            best_ref_values.setdefault(reference, fh.read().strip())
        _, seq = rec
        records_by_reference.setdefault(reference, []).append((sample, seq))

    os.makedirs(args.out_dir, exist_ok=True)
    manifest_path = os.path.join(args.out_dir, 'trees.tsv')
    with open(manifest_path, 'wt') as manifest:
        manifest.write('reference\tdirectory\tsample_count\thas_reference\n')
        for reference in sorted(records_by_reference):
            directory = safe_directory_name(reference)
            virus_dir = os.path.join(args.out_dir, directory)
            records = records_by_reference[reference]

            tree_path = build_tree(records, os.path.join(virus_dir, 'samples_only'), args.mafft_opts)
            print(f'Wrote samples-only tree for {reference} ({len(records)} sample(s)) to {tree_path}')

            ref_seq = find_reference_sequence(reference, best_ref_values[reference], args.reference_dir)
            has_reference = ref_seq is not None
            if has_reference:
                with_ref = [(f'REF_{reference}', ref_seq)] + records
                tree_path = build_tree(with_ref, os.path.join(virus_dir, 'with_reference'), args.mafft_opts)
                print(f'Wrote tree with reference for {reference} to {tree_path}')
            else:
                print(f'Warning: reference FASTA for {reference} not found; skipping tree with reference')
            manifest.write(f'{reference}\t{directory}\t{len(records)}\t{int(has_reference)}\n')

    if not records_by_reference:
        print('No samples with a selected reference were available for phylogeny.')
    print('Wrote tree manifest to', manifest_path)


if __name__ == '__main__':
    main()
