#!/usr/bin/env python3
import argparse
import os

import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
from jinja2 import Environment, FileSystemLoader
import pysam


CIGAR_LABELS = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X",
}


def coverage_bins_and_cigar(bam_path, bin_size=200):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = bam.references
        if not refs:
            return [], [], {}
        ref = refs[0]
        ref_len = bam.get_reference_length(ref)
        bins = (ref_len // bin_size) + 1
        coverage = np.zeros(bins, dtype=np.int64)
        cigar_counts = {label: 0 for label in CIGAR_LABELS.values()}

        for read in bam.fetch(ref):
            if read.is_unmapped:
                continue

            start = max(read.reference_start, 0)
            end = max(read.reference_end or (read.reference_start + read.query_length), start + 1)
            b1 = start // bin_size
            b2 = (end - 1) // bin_size
            coverage[b1:b2 + 1] += 1

            if read.cigartuples:
                for op, length in read.cigartuples:
                    label = CIGAR_LABELS.get(op)
                    if label is not None:
                        cigar_counts[label] += length

    x = [i * bin_size for i in range(len(coverage))]
    return x, coverage.tolist(), cigar_counts


def sample_panels(folder):
    panels = []
    assembly_dir = os.path.join(folder, "assembly")

    for fn in sorted(os.listdir(assembly_dir)):
        if not fn.endswith(".bam"):
            continue
        sample = fn.replace(".bam", "")
        classic_bam = os.path.join(folder, "assembly", f"{sample}.bam")
        snippy_bam = os.path.join(folder, "snippy", sample, "snps.bam")
        if not os.path.exists(snippy_bam):
            continue

        x_c, y_c, cigar_c = coverage_bins_and_cigar(classic_bam)
        x_s, y_s, cigar_s = coverage_bins_and_cigar(snippy_bam)

        if not x_c or not x_s:
            continue

        fig = go.Figure()
        fig.add_trace(go.Bar(name="Classic", x=x_c, y=y_c, marker_color="#5DADE2"))
        fig.add_trace(go.Bar(name="Snippy-style", x=x_s, y=y_s, marker_color="#F1948A"))
        fig.update_layout(
            barmode="overlay",
            title=f"{sample} - read placement blocks on reference",
            xaxis_title="Reference position (bp, bin start)",
            yaxis_title="Coverage per bin",
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=360,
        )

        panels.append(
            {
                "sample": sample,
                "plot_div": plot(fig, output_type="div", include_plotlyjs=False),
                "classic_cigar": cigar_c,
                "snippy_cigar": cigar_s,
            }
        )

    return panels


def main():
    parser = argparse.ArgumentParser(description="Nucleotide placement report")
    parser.add_argument("--folder", required=True)
    parser.add_argument("--install_path", required=True)
    args = parser.parse_args()

    panels = sample_panels(args.folder)

    env = Environment(
        loader=FileSystemLoader(os.path.join(args.install_path, "html_templates")),
        autoescape=False,
    )
    template = env.get_template("nucleotide_placement.html")
    html = template.render(panels=panels)

    out_html = os.path.join(args.folder, "nucleotide_placement_report.html")
    with open(out_html, "wt") as out:
        out.write(html)

    print(f"Wrote nucleotide placement report to {out_html}")


if __name__ == "__main__":
    main()
