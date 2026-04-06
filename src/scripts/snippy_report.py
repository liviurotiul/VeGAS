#!/usr/bin/env python3
import argparse
import os
from itertools import combinations
from collections import defaultdict

import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
from jinja2 import Environment, FileSystemLoader


def parse_vcf_variant_count(vcf_path):
    count = 0
    with open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.strip():
                count += 1
    return count


def parse_fasta_records(fasta_path):
    records = {}
    current_name = None
    current_seq = []

    with open(fasta_path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_name is not None:
        records[current_name] = "".join(current_seq)

    return records


def snp_distance(seq_a, seq_b):
    valid = {"A", "C", "G", "T"}
    dist = 0
    for base_a, base_b in zip(seq_a, seq_b):
        if base_a in valid and base_b in valid and base_a != base_b:
            dist += 1
    return dist


def compute_distance_matrix(records):
    sample_names = list(records.keys())
    n = len(sample_names)
    matrix = np.zeros((n, n), dtype=int)

    for i, j in combinations(range(n), 2):
        d = snp_distance(records[sample_names[i]], records[sample_names[j]])
        matrix[i, j] = d
        matrix[j, i] = d

    return sample_names, matrix


def build_barplot(sample_counts):
    samples = [sample for sample, _ in sample_counts]
    counts = [count for _, count in sample_counts]

    fig = go.Figure(
        data=[go.Bar(x=samples, y=counts, marker_color="#8FB9E3")]
    )
    fig.update_layout(
        title="SNP count per sample",
        xaxis_title="Sample",
        yaxis_title="SNP count",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return plot(fig, output_type="div", include_plotlyjs=False)


def build_heatmap(sample_names, matrix):
    fig = go.Figure(
        data=[
            go.Heatmap(
                z=matrix,
                x=sample_names,
                y=sample_names,
                colorscale="Blues",
                colorbar=dict(title="SNP distance"),
            )
        ]
    )
    fig.update_layout(
        title="Pairwise SNP distance matrix (from core.aln)",
        xaxis_title="Sample",
        yaxis_title="Sample",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return plot(fig, output_type="div", include_plotlyjs=False)


def parse_best_reference(path):
    with open(path, "rt") as handle:
        value = handle.read().strip()
    return os.path.basename(value)


def parse_snps_tab(tab_path, sample_name, reference_name):
    rows = []
    if not os.path.exists(tab_path):
        return rows

    with open(tab_path, "rt") as handle:
        header = None
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            record = dict(zip(header, parts))
            record["sample"] = sample_name
            record["reference"] = reference_name
            rows.append(record)
    return rows


def build_reference_barplot(reference_name, sample_counts):
    samples = [sample for sample, _ in sample_counts]
    counts = [count for _, count in sample_counts]

    fig = go.Figure(
        data=[go.Bar(x=samples, y=counts, marker_color="#5AA469")]
    )
    fig.update_layout(
        title=f"{reference_name} - mutation count per sample",
        xaxis_title="Sample",
        yaxis_title="Mutation count",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return plot(fig, output_type="div", include_plotlyjs=False)


def safe_reference_filename(reference_name):
    return "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in reference_name)


def write_reference_mutation_tsv(path, rows):
    columns = [
        "reference", "sample", "CHROM", "POS", "TYPE", "REF", "ALT", "EVIDENCE",
        "FTYPE", "STRAND", "NT_POS", "AA_POS", "EFFECT", "LOCUS_TAG", "GENE", "PRODUCT",
    ]
    with open(path, "wt") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(col, "")) for col in columns) + "\n")


def build_reference_groups(folder):
    groups = defaultdict(lambda: {"samples": [], "rows": []})
    assembly_dir = os.path.join(folder, "assembly")
    snippy_dir = os.path.join(folder, "snippy")

    if not (os.path.isdir(assembly_dir) and os.path.isdir(snippy_dir)):
        return groups

    for file_name in sorted(os.listdir(assembly_dir)):
        if not file_name.endswith(".best_ref.txt"):
            continue

        sample_name = file_name.replace(".best_ref.txt", "")
        ref_name = parse_best_reference(os.path.join(assembly_dir, file_name))
        tab_path = os.path.join(snippy_dir, sample_name, "snps.tab")
        rows = parse_snps_tab(tab_path, sample_name, ref_name)

        groups[ref_name]["samples"].append(sample_name)
        groups[ref_name]["rows"].extend(rows)

    return groups


def render_reference_reports(folder, install_path):
    groups = build_reference_groups(folder)
    out_dir = os.path.join(folder, "snippy_by_reference")
    os.makedirs(out_dir, exist_ok=True)

    env = Environment(loader=FileSystemLoader(os.path.join(install_path, "html_templates")), autoescape=False)
    template = env.get_template("snippy_by_reference.html")

    index_rows = []
    for reference_name in sorted(groups.keys()):
        sample_names = sorted(groups[reference_name]["samples"])
        mutation_rows = sorted(groups[reference_name]["rows"], key=lambda row: (row.get("sample", ""), int(row.get("POS", 0))))
        sample_counts = []
        for sample_name in sample_names:
            count = sum(1 for row in mutation_rows if row.get("sample") == sample_name)
            sample_counts.append((sample_name, count))

        plot_div = build_reference_barplot(reference_name, sample_counts) if sample_counts else "<p>No mutation rows found.</p>"

        safe_name = safe_reference_filename(reference_name)
        html_name = f"{safe_name}.html"
        tsv_name = f"{safe_name}.tsv"
        html_path = os.path.join(out_dir, html_name)
        tsv_path = os.path.join(out_dir, tsv_name)

        write_reference_mutation_tsv(tsv_path, mutation_rows)

        rendered = template.render(
            reference_name=reference_name,
            sample_names=sample_names,
            sample_counts=sample_counts,
            mutation_rows=mutation_rows,
            plot_div=plot_div,
            tsv_name=tsv_name,
        )
        with open(html_path, "wt") as handle:
            handle.write(rendered)

        index_rows.append(
            {
                "reference_name": reference_name,
                "sample_count": len(sample_names),
                "mutation_count": len(mutation_rows),
                "html_name": html_name,
                "tsv_name": tsv_name,
            }
        )

    index_path = os.path.join(out_dir, "index.html")
    with open(index_path, "wt") as handle:
        handle.write(
            "<!DOCTYPE html><html><head><meta charset='UTF-8'><title>Snippy By Reference</title>"
            "<link rel='stylesheet' href='https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css'>"
            "</head><body><div class='container'><h1 class='mt-4'>Snippy reports grouped by reference</h1>"
            "<p>Each report groups samples that share the same selected reference and lists mutations per sample.</p>"
            "<table class='table table-striped table-sm'><thead><tr><th>Reference</th><th>Samples</th><th>Mutations</th><th>HTML</th><th>TSV</th></tr></thead><tbody>"
        )
        for row in index_rows:
            handle.write(
                f"<tr><td>{row['reference_name']}</td><td>{row['sample_count']}</td><td>{row['mutation_count']}</td>"
                f"<td><a href='{row['html_name']}'>open</a></td><td><a href='{row['tsv_name']}'>download</a></td></tr>"
            )
        handle.write("</tbody></table></div></body></html>")

    return out_dir, index_rows


def main():
    parser = argparse.ArgumentParser(description="Generate Snippy summary report")
    parser.add_argument("--folder", required=True, help="Working directory")
    parser.add_argument("--install_path", required=True, help="Install/src path")
    args = parser.parse_args()

    snippy_dir = os.path.join(args.folder, "snippy")
    core_aln = os.path.join(snippy_dir, "core.aln")

    sample_counts = []
    if os.path.isdir(snippy_dir):
        for name in sorted(os.listdir(snippy_dir)):
            sample_dir = os.path.join(snippy_dir, name)
            if not os.path.isdir(sample_dir):
                continue
            vcf_path = os.path.join(sample_dir, "snps.vcf")
            if os.path.exists(vcf_path):
                sample_counts.append((name, parse_vcf_variant_count(vcf_path)))

    sample_names = []
    matrix = np.zeros((0, 0), dtype=int)
    if os.path.exists(core_aln):
        records = parse_fasta_records(core_aln)
        if records:
            sample_names, matrix = compute_distance_matrix(records)

    barplot_div = build_barplot(sample_counts) if sample_counts else "<p>No per-sample snps.vcf files found.</p>"
    heatmap_div = (
        build_heatmap(sample_names, matrix)
        if len(sample_names) > 0
        else "<p>No core.aln found or alignment is empty.</p>"
    )

    reference_report_dir, reference_rows = render_reference_reports(args.folder, args.install_path)

    env = Environment(loader=FileSystemLoader(os.path.join(args.install_path, "html_templates")), autoescape=False)
    template = env.get_template("snippy.html")

    rendered = template.render(
        sample_counts=sample_counts,
        barplot_div=barplot_div,
        heatmap_div=heatmap_div,
        reference_rows=reference_rows,
        reference_report_dir=os.path.basename(reference_report_dir),
    )

    out_path = os.path.join(args.folder, "snippy_report.html")
    with open(out_path, "wt") as handle:
        handle.write(rendered)

    print(f"Wrote Snippy report to {out_path}")
    print(f"Wrote grouped-by-reference reports to {reference_report_dir}")


if __name__ == "__main__":
    main()
