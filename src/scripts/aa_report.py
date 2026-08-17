#!/usr/bin/env python3
"""
aa_report.py — Compară aminoacizii: Referință vs Snippy vs Clasic.

Generează aa_report.html în <folder>/.

Utilizare:
  python aa_report.py --folder data/output --reference_dir data/reference
"""
import argparse
import os
import sys
from glob import glob
from html import escape

# Asigură că directorul scriptului este în Python path (import local pathogen_db)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)
from pathogen_db import PATHOGEN_DB, detect_pathogen

# ── Numele complete ale aminoacizilor ──────────────────────────────────────
AA_NAMES = {
    "A": "Alanine (Ala)",
    "R": "Arginine (Arg)",
    "N": "Asparagine (Asn)",
    "D": "Aspartic acid (Asp)",
    "C": "Cysteine (Cys)",
    "Q": "Glutamine (Gln)",
    "E": "Glutamic acid (Glu)",
    "G": "Glycine (Gly)",
    "H": "Histidine (His)",
    "I": "Isoleucine (Ile)",
    "L": "Leucine (Leu)",
    "K": "Lysine (Lys)",
    "M": "Methionine (Met)",
    "F": "Phenylalanine (Phe)",
    "P": "Proline (Pro)",
    "S": "Serine (Ser)",
    "T": "Threonine (Thr)",
    "W": "Tryptophan (Trp)",
    "Y": "Tyrosine (Tyr)",
    "V": "Valine (Val)",
    "*": "Stop codon",
    "?": "Nedeterminat (lipsă acoperire)",
    "X": "Necunoscut",
}

# ── Codul genetic standard ──────────────────────────────────────────────────
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate(nt):
    """Returnează șirul de aminoacizi."""
    nt = nt.upper()
    aa = []
    for i in range(0, len(nt) - 2, 3):
        codon = nt[i:i+3]
        if "N" in codon or "-" in codon:
            aa.append("?")
        else:
            aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)

def get_codons(nt):
    """Returnează lista de codoni (câte 3 NT) pentru fiecare poziție AA."""
    nt = nt.upper()
    codons = []
    for i in range(0, len(nt) - 2, 3):
        codon = nt[i:i+3]
        codons.append(codon)
    return codons

def parse_fasta(path):
    records = {}
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name:
        records[name] = "".join(parts)
    return records

def parse_gff_cds(path):
    features = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            chrom, _, _, start, end, _, strand, phase_str, attrs = cols[:9]
            info = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k.strip()] = v.strip()
            features.append({
                "name":    info.get("Name", info.get("ID", "?")),
                "product": info.get("product", ""),
                "chrom":   chrom,
                "start":   int(start),   # 1-based
                "end":     int(end),     # 1-based inclusiv
                "strand":  strand,
                "phase":   int(phase_str) if phase_str in ('0', '1', '2') else 0,
            })
    return features

def extract_cds(genome_dict, feat):
    """Extrage secvența CDS din dicționarul genomic."""
    seq = genome_dict.get(feat["chrom"], "")
    if not seq:
        # încearcă primul contig disponibil (pentru assembly clasic)
        if genome_dict:
            seq = next(iter(genome_dict.values()))
    if not seq:
        return ""
    nt = seq[feat["start"] - 1 : feat["end"]]
    if feat["strand"] == "-":
        comp = str.maketrans("ACGTN", "TGCAN")
        nt = nt.translate(comp)[::-1]
    phase = feat.get("phase", 0)
    if phase:
        nt = nt[phase:]
    return nt

def find_samples(folder):
    snippy_dir = os.path.join(folder, "snippy")
    samples = []
    for entry in sorted(os.listdir(snippy_dir)):
        if os.path.isfile(os.path.join(snippy_dir, entry, "snps.consensus.fa")):
            samples.append(entry)
    return samples

def read_sample_best_refs(folder, samples):
    """Citește assembly/{sample}.best_ref.txt și returnează {sample: ref_name}.
    Extrage doar numele de bază (fără prefix ref_indexes/ etc.)."""
    result = {}
    for s in samples:
        best_ref_path = os.path.join(folder, "assembly", f"{s}.best_ref.txt")
        if os.path.isfile(best_ref_path):
            with open(best_ref_path) as fh:
                raw = fh.read().strip()
            result[s] = os.path.basename(raw)
        else:
            result[s] = None
    return result

# ── HTML helpers ────────────────────────────────────────────────────────────
CSS = """
<style>
  body{font-family:Arial,sans-serif;margin:30px;background:#f5f5f5;}
  h1{color:#2c3e50;}
  h2{color:#2c3e50;border-bottom:2px solid #ddd;padding-bottom:4px;margin-top:40px;}
  h3{color:#555;margin-top:28px;}
  .summary{display:flex;gap:14px;margin:18px 0;flex-wrap:wrap;}
  .card{background:#fff;border-radius:8px;padding:12px 20px;
        box-shadow:0 2px 6px rgba(0,0,0,.1);text-align:center;min-width:100px;}
  .card .num{font-size:1.8em;font-weight:bold;}
  .card .lbl{font-size:.8em;color:#888;}
  .tbl{border-collapse:collapse;font-size:.87em;background:#fff;
       box-shadow:0 2px 6px rgba(0,0,0,.1);margin-bottom:8px;}
  .tbl th{background:#2c3e50;color:#fff;padding:7px 10px;
          text-align:center;position:sticky;top:0;z-index:1;}
  .tbl td{padding:5px 9px;border:1px solid #e5e5e5;font-family:monospace;}
  .tbl td.pos{background:#f0f4f8;color:#555;font-weight:bold;text-align:right;}
  .tbl td.ref{background:#eaf4fb;color:#2c3e50;}
  .tbl td.same{color:#ccc;}
  .tbl td.diff{background:#e74c3c;color:#fff;font-weight:bold;}
  .tbl td.stop{background:#8e44ad;color:#fff;font-weight:bold;}
  .tbl td.miss{background:#ffeaa7;color:#333;}
  .tbl td.unk{color:#aaa;font-style:italic;}
  .tbl tr:hover td{filter:brightness(.93);}
  .legend{display:flex;gap:8px;flex-wrap:wrap;align-items:center;
          background:#fff;border-radius:8px;padding:10px 16px;
          box-shadow:0 2px 6px rgba(0,0,0,.1);margin-bottom:20px;font-size:.85em;}
  .leg{padding:3px 10px;border-radius:4px;}
  .info{background:#eaf4fb;border-left:4px solid #3498db;
        padding:10px 16px;border-radius:4px;margin-bottom:20px;
        font-size:.9em;color:#444;}
  .nav{background:#fff;border-radius:8px;padding:10px 16px;
       box-shadow:0 2px 6px rgba(0,0,0,.1);margin-bottom:20px;}
  details summary{cursor:pointer;font-weight:bold;color:#2c3e50;}
  details{margin-bottom:8px;}
  .wt{color:#27ae60;margin:8px 0;}
  .overflow{overflow-x:auto;}
  /* tabel cod genetic */
  .aatbl{border-collapse:collapse;font-size:.86em;background:#fff;
         box-shadow:0 2px 6px rgba(0,0,0,.1);}
  .aatbl th{background:#34495e;color:#fff;padding:6px 12px;text-align:left;}
  .aatbl td{padding:5px 12px;border:1px solid #e5e5e5;}
  .aatbl td.letter{font-family:monospace;font-size:1.1em;font-weight:bold;
                    text-align:center;background:#f8f8f8;}
  .aatbl td.codons{font-family:monospace;color:#555;font-size:.9em;}
</style>
"""

# ── Tabel legendă aminoacizi ─────────────────────────────────────────────────
# Numele complet + codoni pentru fiecare AA (grupați frumos)
AA_DETAILS = [
    ("A","Alanine (Ala)",          "GCT GCC GCA GCG"),
    ("R","Arginine (Arg)",          "CGT CGC CGA CGG AGA AGG"),
    ("N","Asparagine (Asn)",        "AAT AAC"),
    ("D","Aspartic acid (Asp)",     "GAT GAC"),
    ("C","Cysteine (Cys)",          "TGT TGC"),
    ("Q","Glutamine (Gln)",         "CAA CAG"),
    ("E","Glutamic acid (Glu)",     "GAA GAG"),
    ("G","Glycine (Gly)",           "GGT GGC GGA GGG"),
    ("H","Histidine (His)",         "CAT CAC"),
    ("I","Isoleucine (Ile)",        "ATT ATC ATA"),
    ("L","Leucine (Leu)",           "TTA TTG CTT CTC CTA CTG"),
    ("K","Lysine (Lys)",            "AAA AAG"),
    ("M","Methionine (Met) / START","ATG"),
    ("F","Phenylalanine (Phe)",     "TTT TTC"),
    ("P","Proline (Pro)",           "CCT CCC CCA CCG"),
    ("S","Serine (Ser)",            "TCT TCC TCA TCG AGT AGC"),
    ("T","Threonine (Thr)",         "ACT ACC ACA ACG"),
    ("W","Tryptophan (Trp)",        "TGG"),
    ("Y","Tyrosine (Tyr)",          "TAT TAC"),
    ("V","Valine (Val)",            "GTT GTC GTA GTG"),
    ("*","Stop codon",              "TAA TAG TGA"),
]

def aa_legend_html():
    rows = ""
    for letter, name, codons in AA_DETAILS:
        rows += (
            f"<tr>"
            f'<td class="letter">{escape(letter)}</td>'
            f"<td>{escape(name)}</td>"
            f'<td class="codons">{escape(codons)}</td>'
            f"</tr>\n"
        )
    return (
        '<div style="overflow-x:auto">'
        '<table class="aatbl">'
        "<thead><tr><th>Cod</th><th>Aminoacid</th><th>Codoni (triplete NT)</th></tr></thead>"
        f"<tbody>{rows}</tbody>"
        "</table></div>"
    )

# ═══════════════ Rezistență la antivirale / antifungice ══════════════════════

def check_resistance(pathogen_id, samples, sample_snippy, sample_classic, features):
    """Detectează mutații de rezistență pentru patogenul dat.

    Returnează {sample: {gene_canonical: [hit_dict, ...]}}.
    hit_dict are cheile: pos, wt, mut, level, drug_class, drugs, note, gene, method.
    """
    if pathogen_id not in PATHOGEN_DB:
        return {s: {} for s in samples}

    pdata     = PATHOGEN_DB[pathogen_id]
    aliases   = pdata["gene_aliases"]
    mutations = pdata["mutations"]

    # Mapează feature-urile GFF la gene canonice
    gene_features = {}
    for feat in features:
        canonical = aliases.get(feat["name"].upper())
        if canonical:
            gene_features[canonical] = feat

    results = {}
    for s in samples:
        results[s] = {}
        for canon, feat in gene_features.items():
            snp_nt  = extract_cds(sample_snippy[s], feat)
            cls_nt  = extract_cds(sample_classic.get(s, {}), feat)
            snp_seq = translate(snp_nt)
            cls_seq = translate(cls_nt) if cls_nt else ""

            hits = []
            for (gene, pos, wt, alts_str, drug_class, drugs, level, note) in mutations:
                if gene != canon:
                    continue
                idx      = pos - 1
                alt_list = [a.strip() for a in alts_str.split("/")]
                for method, seq in [("Snippy", snp_seq), ("Clasic", cls_seq)]:
                    if idx >= len(seq) or seq[idx] in ("?", "X"):
                        continue
                    if seq[idx] in alt_list:
                        hits.append({
                            "pos": pos, "wt": wt, "mut": seq[idx],
                            "level": level, "drug_class": drug_class,
                            "drugs": drugs, "note": note,
                            "gene": canon, "method": method,
                        })
            if hits:
                results[s][canon] = hits
    return results


def build_resistance_html(pathogen_id, resistance_results, reason, confidence):
    """Construiește secțiunea HTML de rezistență pentru orice patogen."""
    if pathogen_id is None or pathogen_id not in PATHOGEN_DB:
        # Nu afișăm nicio secțiune de rezistență dacă referința nu e un patogen cunoscut
        # (ex: hMPV, SARS-CoV-2 fără bază de date, etc.)
        return ""

    pdata       = PATHOGEN_DB[pathogen_id]
    drug_classes = pdata["drug_classes"]
    dc_order    = pdata["drug_class_order"]
    total       = sum(len(h) for sd in resistance_results.values() for h in sd.values())

    conf_html = (
        f'<span style="background:#27ae60;color:#fff;padding:2px 10px;'
        f'border-radius:12px;font-size:.82em">Confirmat: {escape(pdata["short_name"])}</span>'
        if confidence == "high" else
        f'<span style="background:#f39c12;color:#fff;padding:2px 10px;'
        f'border-radius:12px;font-size:.82em">Probabil: {escape(pdata["short_name"])}</span>'
    )

    blocks = []
    for sample in sorted(resistance_results):
        gene_data = resistance_results[sample]
        short     = (sample[:22] + "…") if len(sample) > 24 else sample

        if not gene_data:
            blocks.append(
                f'<div style="background:#fff;border-radius:8px;padding:14px 18px;'
                f'box-shadow:0 2px 6px rgba(0,0,0,.1);margin-bottom:14px">'
                f'<b>{escape(short)}</b> — '
                f'<span style="color:#27ae60">✔ Nicio mutație de rezistență detectată</span>'
                f'</div>'
            )
            continue

        by_class = {}
        for gene, hits in gene_data.items():
            for h in hits:
                by_class.setdefault(h["drug_class"], []).append(h)

        class_blocks = []
        for dc in dc_order:
            if dc not in by_class:
                continue
            dc_info  = drug_classes[dc]
            bg, fg   = dc_info["bg"], dc_info["fg"]
            rows = ""
            for h in sorted(by_class[dc], key=lambda x: (x["gene"], x["pos"])):
                notation = f'{h["wt"]}{h["pos"]}{h["mut"]}'
                if h["level"] == "major":
                    level_html = ('<span style="background:#e74c3c;color:#fff;padding:1px 7px;'
                                  'border-radius:10px;font-size:.78em">Major</span>')
                else:
                    level_html = ('<span style="background:#f39c12;color:#fff;padding:1px 7px;'
                                  'border-radius:10px;font-size:.78em">Minor/Moderate</span>')
                note_html = (
                    f' <span style="color:#aaa;font-size:.83em">({escape(h["note"])})</span>'
                    if h["note"] else ""
                )
                rows += (
                    f"<tr>"
                    f'<td style="font-family:monospace;font-weight:bold;padding:4px 8px">{escape(h["gene"])}</td>'
                    f'<td style="text-align:right;padding:4px 8px">{h["pos"]}</td>'
                    f'<td style="font-family:monospace;color:#888;padding:4px 8px">{escape(h["wt"])}</td>'
                    f'<td style="font-family:monospace;font-weight:bold;color:#e74c3c;padding:4px 8px">{escape(h["mut"])}</td>'
                    f'<td style="font-family:monospace;padding:4px 8px">{escape(notation)}</td>'
                    f'<td style="padding:4px 8px">{level_html}</td>'
                    f'<td style="font-size:.85em;padding:4px 8px">{escape(", ".join(h["drugs"]))}{note_html}</td>'
                    f'<td style="font-size:.85em;color:#555;padding:4px 8px">{escape(h["method"])}</td>'
                    f"</tr>"
                )
            class_blocks.append(
                f'<div style="margin-bottom:12px">'
                f'<div style="background:{bg};color:{fg};padding:5px 14px;'
                f'border-radius:6px 6px 0 0;font-weight:bold;font-size:.88em">'
                f'{escape(dc_info["name"])}</div>'
                f'<table style="width:100%;border-collapse:collapse;border:1px solid #ddd;'
                f'font-size:.86em;background:#fff">'
                f'<thead style="background:#f5f5f5"><tr>'
                f'<th style="padding:5px 8px">Genă</th><th style="padding:5px 8px">Pos</th>'
                f'<th style="padding:5px 8px">WT</th><th style="padding:5px 8px">Mut</th>'
                f'<th style="padding:5px 8px">Notație</th><th style="padding:5px 8px">Nivel</th>'
                f'<th style="padding:5px 8px">Droguri afectate</th>'
                f'<th style="padding:5px 8px">Metodă</th>'
                f'</tr></thead><tbody>{rows}</tbody></table></div>'
            )

        blocks.append(
            f'<div style="background:#fff;border-radius:8px;padding:16px 20px;'
            f'box-shadow:0 2px 6px rgba(0,0,0,.1);margin-bottom:18px">'
            f'<h3 style="margin-top:0;color:#2c3e50">{escape(short)}</h3>'
            + "".join(class_blocks)
            + "</div>"
        )

    return (
        '<hr style="margin-top:50px">'
        f'<h2 id="resistance">Analiză rezistență — {escape(pdata["display_name"])}</h2>'
        f'<div class="info">{conf_html} &nbsp; {escape(reason)}<br><br>'
        f'Mutații comparate cu <b>{escape(pdata["source"])}</b>. '
        'Se raportează mutațiile față de wild-type în secvențele Snippy și Clasic.<br>'
        f'<a href="{escape(pdata["source_url"])}" target="_blank" style="color:#2980b9">'
        f'{escape(pdata["source_url"])}</a>'
        '</div>'
        f'<div class="summary">'
        f'<div class="card"><div class="num" style="color:#e74c3c">{total}</div>'
        f'<div class="lbl">Mutații rezistență găsite</div></div>'
        f'<div class="card"><div class="num">{len(resistance_results)}</div>'
        f'<div class="lbl">Probe analizate</div></div></div>'
        + "".join(blocks)
    )


def build_resistance_page(pathogen_id, resistance_results, reason, confidence):
    """Construiește o pagină HTML completă standalone pentru raportul de rezistență HIV."""
    if pathogen_id is None or pathogen_id not in PATHOGEN_DB:
        return None

    pdata        = PATHOGEN_DB[pathogen_id]
    drug_classes = pdata["drug_classes"]
    dc_order     = pdata["drug_class_order"]
    total        = sum(len(h) for sd in resistance_results.values() for h in sd.values())

    # Rezumat global per clasă de drog
    global_by_class = {}
    for sd in resistance_results.values():
        for gene, hits in sd.items():
            for h in hits:
                global_by_class.setdefault(h["drug_class"], set()).add(
                    f'{h["wt"]}{h["pos"]}{h["mut"]}'
                )

    summary_cards = ""
    summary_cards += (
        f'<div class="card"><div class="num" style="color:#e74c3c">{total}</div>'
        f'<div class="lbl">Mutații rezistență găsite</div></div>'
        f'<div class="card"><div class="num">{len(resistance_results)}</div>'
        f'<div class="lbl">Probe analizate</div></div>'
    )
    for dc in dc_order:
        n = len(global_by_class.get(dc, []))
        if n:
            bg = pdata["drug_classes"][dc]["bg"]
            summary_cards += (
                f'<div class="card" style="border-top:4px solid {bg}">'
                f'<div class="num" style="color:{bg}">{n}</div>'
                f'<div class="lbl">{escape(pdata["drug_classes"][dc]["name"])}</div></div>'
            )

    # Secțiune per probă
    sample_sections = []
    for sample in sorted(resistance_results):
        gene_data = resistance_results[sample]
        short = (sample[:40] + "…") if len(sample) > 42 else sample

        if not gene_data:
            sample_sections.append(
                f'<div class="sample-block">'
                f'<h2>{escape(short)}</h2>'
                f'<p style="color:#27ae60;font-size:1.05em">'
                f'✔ Nicio mutație de rezistență detectată față de wild-type.</p>'
                f'</div>'
            )
            continue

        by_class = {}
        for gene, hits in gene_data.items():
            for h in hits:
                by_class.setdefault(h["drug_class"], []).append(h)

        class_blocks = []
        for dc in dc_order:
            if dc not in by_class:
                continue
            dc_info = drug_classes[dc]
            bg, fg  = dc_info["bg"], dc_info["fg"]
            rows = ""
            for h in sorted(by_class[dc], key=lambda x: (x["gene"], x["pos"])):
                notation = f'{h["wt"]}{h["pos"]}{h["mut"]}'
                level_style = (
                    "background:#e74c3c;color:#fff" if h["level"] == "major"
                    else "background:#f39c12;color:#fff"
                )
                level_label = "Major" if h["level"] == "major" else "Minor/Moderate"
                note_html = (
                    f'<br><span style="color:#888;font-size:.82em">{escape(h["note"])}</span>'
                    if h["note"] else ""
                )
                rows += (
                    f"<tr>"
                    f'<td class="mono bold">{escape(h["gene"])}</td>'
                    f'<td class="right">{h["pos"]}</td>'
                    f'<td class="mono wt">{escape(h["wt"])}</td>'
                    f'<td class="mono mut">{escape(h["mut"])}</td>'
                    f'<td class="mono notation">{escape(notation)}</td>'
                    f'<td><span class="badge" style="{level_style}">{level_label}</span></td>'
                    f'<td class="drugs">{escape(", ".join(h["drugs"]))}{note_html}</td>'
                    f'<td class="method">{escape(h["method"])}</td>'
                    f"</tr>"
                )
            class_blocks.append(
                f'<div class="drug-class-block">'
                f'<div class="dc-header" style="background:{bg};color:{fg}">'
                f'{escape(dc_info["name"])}</div>'
                f'<table class="rtbl"><thead><tr>'
                f'<th>Genă</th><th>Pos</th><th>WT</th><th>Mut</th>'
                f'<th>Notație</th><th>Nivel</th><th>Droguri afectate</th><th>Metodă</th>'
                f'</tr></thead><tbody>{rows}</tbody></table>'
                f'</div>'
            )

        sample_sections.append(
            f'<div class="sample-block">'
            f'<h2>{escape(short)}</h2>'
            + "".join(class_blocks) +
            f'</div>'
        )

    conf_html = (
        f'<span class="badge" style="background:#27ae60;color:#fff">'
        f'Confirmat: {escape(pdata["short_name"])}</span>'
        if confidence == "high" else
        f'<span class="badge" style="background:#f39c12;color:#fff">'
        f'Probabil: {escape(pdata["short_name"])}</span>'
    )

    body = f"""
<h1>VeGAS — Raport Rezistență HIV-1</h1>
<div class="info-box">
  {conf_html} &nbsp; {escape(reason)}<br><br>
  Mutații comparate cu <b>{escape(pdata["source"])}</b>.<br>
  Sursă oficială:
  <a href="{escape(pdata["source_url"])}" target="_blank">{escape(pdata["source_url"])}</a><br><br>
  <b>Metodologie:</b> Se identifică mutațiile față de wild-type în aminoacizii genelor
  PR, RT și IN reconstruite prin Snippy și prin aliniere clasică Bowtie2+bcftools.
  Sunt raportate <b>doar mutațiile care apar în baza de date Stanford HIVdb</b> cu relevanță
  clinică dovedită pentru rezistența la antiretrovirale.
  <a href="aa_report.html" style="color:#2980b9">← Înapoi la raportul de aminoacizi</a>
</div>

<div class="summary">{summary_cards}</div>

<h2 style="border-bottom:2px solid #ddd;padding-bottom:6px">Legende niveluri de rezistență</h2>
<div class="legend-box">
  <span class="badge" style="background:#e74c3c;color:#fff">Major</span>
  Mutație care conferă rezistență clinică semnificativă — de obicei suficientă singură
  pentru a reduce eficiența tratamentului.
  &nbsp;&nbsp;
  <span class="badge" style="background:#f39c12;color:#fff">Minor/Moderate</span>
  Contribuție parțială sau intermediară — relevantă mai ales în combinație cu alte mutații.
</div>

{"".join(sample_sections)}
"""

    css_page = """
<style>
  body { font-family: Arial, sans-serif; margin: 30px; background: #f5f5f5; }
  h1 { color: #2c3e50; }
  h2 { color: #2c3e50; margin-top: 36px; }
  .info-box { background: #eaf4fb; border-left: 4px solid #3498db;
              padding: 14px 18px; border-radius: 4px; margin-bottom: 24px;
              font-size: .93em; color: #444; }
  .legend-box { background: #fff; border-radius: 8px; padding: 12px 18px;
                box-shadow: 0 2px 6px rgba(0,0,0,.1); margin-bottom: 28px;
                font-size: .9em; }
  .summary { display: flex; gap: 14px; margin: 18px 0; flex-wrap: wrap; }
  .card { background: #fff; border-radius: 8px; padding: 12px 20px;
          box-shadow: 0 2px 6px rgba(0,0,0,.1); text-align: center; min-width: 110px; }
  .card .num { font-size: 1.8em; font-weight: bold; }
  .card .lbl { font-size: .8em; color: #888; }
  .badge { padding: 2px 10px; border-radius: 12px; font-size: .83em;
           display: inline-block; }
  .sample-block { background: #fff; border-radius: 10px; padding: 22px 26px;
                  box-shadow: 0 2px 8px rgba(0,0,0,.1); margin-bottom: 28px; }
  .sample-block h2 { margin-top: 0; font-size: 1.15em; color: #2c3e50; }
  .drug-class-block { margin-bottom: 18px; }
  .dc-header { padding: 6px 14px; border-radius: 6px 6px 0 0;
               font-weight: bold; font-size: .9em; }
  .rtbl { width: 100%; border-collapse: collapse; font-size: .87em;
          border: 1px solid #ddd; }
  .rtbl th { background: #f5f5f5; padding: 6px 10px; text-align: left;
             border-bottom: 2px solid #ddd; }
  .rtbl td { padding: 5px 10px; border: 1px solid #eee; vertical-align: top; }
  .rtbl tr:hover td { background: #fafafa; }
  td.mono { font-family: monospace; }
  td.bold { font-weight: bold; }
  td.right { text-align: right; }
  td.wt { color: #27ae60; }
  td.mut { color: #e74c3c; font-weight: bold; }
  td.notation { font-weight: bold; }
  td.drugs { font-size: .85em; }
  td.method { color: #888; font-size: .85em; }
</style>
"""

    return f"""<!DOCTYPE html>
<html lang="ro">
<head>
<meta charset="UTF-8">
<title>VeGAS — Raport Rezistență HIV-1</title>
{css_page}
</head>
<body>
{body}
</body>
</html>
"""


# ═════════════════════════════════════════════════════════════════════════════

def aa_tooltip(aa, codon=""):
    """Construiește textul tooltip: codon + numele complet al aminoacidului."""
    name = AA_NAMES.get(aa, aa)
    if codon:
        return f"{codon} → {aa} ({name})"
    return f"{aa} ({name})"

def cell(aa, ref_aa, codon=""):
    a = escape(aa)
    tip = escape(aa_tooltip(aa, codon))
    if aa == "?":
        return f'<td class="unk" title="{tip}">{a}</td>'
    if aa == ref_aa:
        return f'<td class="same" title="{tip}">{a}</td>'
    if aa == "*":
        return f'<td class="stop" title="{tip}">{a}</td>'
    return f'<td class="diff" title="{tip}">{a}</td>'

def build_report(samples, ref_genome, sample_snippy, sample_classic, features, ref_name, full_doc=True):
    total_missense = 0
    gene_sections = []
    nav_links = []

    for feat in features:
        gene = feat["name"]
        product = feat["product"]
        anchor = f"g-{gene}"
        nav_links.append((anchor, gene))

        # Traduce referința
        ref_nt = extract_cds(ref_genome, feat)
        ref_aa = translate(ref_nt)
        ref_codons = get_codons(ref_nt)

        # Traduce per probă: snippy și clasic
        snp_aas = {}
        cls_aas = {}
        snp_codons = {}
        cls_codons = {}
        for s in samples:
            snt = extract_cds(sample_snippy[s], feat)
            cnt = extract_cds(sample_classic[s], feat)
            snp_aas[s] = translate(snt)
            cls_aas[s] = translate(cnt)
            snp_codons[s] = get_codons(snt)
            cls_codons[s] = get_codons(cnt)

        # Găsește pozițiile cu vreo diferență față de REF (în oricare probă/metodă)
        max_len = max([len(ref_aa)] +
                      [len(snp_aas[s]) for s in samples] +
                      [len(cls_aas[s]) for s in samples])
        diff_pos = []
        for i in range(max_len):
            r = ref_aa[i] if i < len(ref_aa) else "?"
            for s in samples:
                sn = snp_aas[s][i] if i < len(snp_aas[s]) else "?"
                cl = cls_aas[s][i] if i < len(cls_aas[s]) else "?"
                if (sn != r and sn != "?") or (cl != r and cl != "?"):
                    diff_pos.append(i)
                    total_missense += 1
                    break

        # Header tabel — coloane duble per probă (Snippy | Clasic)
        sample_headers = ""
        for s in samples:
            short = s[:14] + "…" if len(s) > 16 else s
            sample_headers += (
                f'<th colspan="2" title="{escape(s)}">{escape(short)}</th>'
            )
        sub_headers = "".join("<th>Snippy</th><th>Clasic</th>" for _ in samples)
        thead = (
            f"<tr><th>Pos AA</th><th>REF</th>{sample_headers}</tr>"
            f"<tr><th></th><th></th>{sub_headers}</tr>"
        )

        if not diff_pos:
            gene_sections.append(
                f'<div id="{escape(anchor)}">'
                f'<h3>Gen <b>{escape(gene)}</b>'
                f'<span style="font-weight:normal;color:#888;font-size:.85em"> — {escape(product)}</span></h3>'
                f'<p class="wt">✔ Proteina identică cu REF în toate probele și metode.</p>'
                f'</div>'
            )
            continue

        rows_html = []
        for i in diff_pos:
            r   = ref_aa[i]     if i < len(ref_aa)     else "?"
            rc  = ref_codons[i] if i < len(ref_codons)  else ""
            ref_tip = escape(aa_tooltip(r, rc))
            row = (f'<td class="pos">{i+1}</td>'
                   f'<td class="ref" title="{ref_tip}">{escape(r)}</td>')
            for s in samples:
                sn = snp_aas[s][i]    if i < len(snp_aas[s])    else "?"
                sc = snp_codons[s][i] if i < len(snp_codons[s]) else ""
                cl = cls_aas[s][i]    if i < len(cls_aas[s])    else "?"
                cc = cls_codons[s][i] if i < len(cls_codons[s]) else ""
                row += cell(sn, r, sc) + cell(cl, r, cc)
            rows_html.append(f"<tr>{row}</tr>")

        gene_sections.append(
            f'<div id="{escape(anchor)}">'
            f'<h3>Gen <b>{escape(gene)}</b>'
            f'<span style="font-weight:normal;color:#888;font-size:.85em"> — {escape(product)}'
            f' &nbsp;|&nbsp; {len(diff_pos)} poziții diferite / {len(ref_aa)} AA total</span></h3>'
            f'<div class="overflow">'
            f'<table class="tbl"><thead>{thead}</thead><tbody>'
            + "\n".join(rows_html) +
            f'</tbody></table></div></div>'
        )

    nav_html = " ".join(
        f'<a href="#{a}" style="color:#2c3e50;margin-right:12px">{escape(g)}</a>'
        for a, g in nav_links
    )
    n_genes = len(features)
    n_samples = len(samples)

    body = f"""
<h1>VeGAS — Aminoacizi: REF vs Snippy vs Clasic</h1>
<p>Referință: <b>{escape(ref_name)}</b> &nbsp;|&nbsp; {n_samples} probe &nbsp;|&nbsp; {n_genes} gene CDS</p>

<div class="summary">
  <div class="card"><div class="num" style="color:#e74c3c">{total_missense}</div>
    <div class="lbl">Poziții AA diferite față de REF</div></div>
  <div class="card"><div class="num">{n_genes}</div><div class="lbl">Gene CDS</div></div>
  <div class="card"><div class="num">{n_samples}</div><div class="lbl">Probe</div></div>
</div>

<div class="info">
  <b>Metodologie:</b> Fiecare regiune CDS este extrasă din cele 3 surse și tradusă cu codul genetic standard.
  <b>REF</b> = referința originală &nbsp;|&nbsp;
  <b>Snippy</b> = consens reconstruit de Snippy (<code>snps.consensus.fa</code>) &nbsp;|&nbsp;
  <b>Clasic</b> = assembly prin aliniere Bowtie2+bcftools (<code>assembly/{"{sample}"}.fasta</code>).
  Se afișează <b>doar pozițiile unde cel puțin o probă diferă față de REF</b>. Poziția <code>?</code> = lipsă acoperire (N în secvență).
</div>

<div class="legend">
  <b>Legendă:</b>
  <span class="leg" style="background:#eaf4fb;color:#2c3e50">REF</span>
  <span class="leg" style="background:#e74c3c;color:#fff">Missense / diferit față de REF</span>
  <span class="leg" style="background:#8e44ad;color:#fff">Stop câștigat (*)</span>
  <span class="leg" style="color:#ccc;border:1px solid #eee">Identic cu REF</span>
  <span class="leg" style="color:#aaa;font-style:italic">? lipsă acoperire</span>
</div>

<div class="nav"><b>Navigare:</b> {nav_html}</div>

{"".join(gene_sections)}
<hr style="margin-top:50px">
<h2 id="aa-table">Tabel aminoacizi — codul genetic standard</h2>
{aa_legend_html()}"""

    if not full_doc:
        return body
    return f"""<!DOCTYPE html>
<html lang="ro">
<head>
<meta charset="UTF-8">
<title>VeGAS — Aminoacizi REF vs Snippy vs Clasic</title>
{CSS}
</head>
<body>
{body}
</body>
</html>
"""

# ── Main ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", required=True,
                        help="Director output pipeline (conține snippy/, assembly/).")
    parser.add_argument("--install_path", default=None,
                        help="Neutilizat, păstrat pentru compatibilitate cu Snakefile.")
    parser.add_argument("--reference_dir", default=None,
                        help="Director cu .fasta + .gff de referință.")
    parser.add_argument("--force-hiv", action="store_true", default=False,
                        help="Forțează analiza HIVdb indiferent de virus (demo/test).")
    args = parser.parse_args()

    folder = args.folder
    ref_dir = args.reference_dir or os.path.join(
        os.path.dirname(os.path.abspath(folder)), "reference"
    )

    # Construiește index: ref_name -> (fasta_path, gff_path) pentru toate .gff disponibile
    gff_files = glob(os.path.join(ref_dir, "*.gff"))
    if not gff_files:
        print(f"[aa_report] EROARE: niciun .gff în {ref_dir}")
        raise SystemExit(1)

    ref_index = {}
    for gff_path in gff_files:
        ref_name = os.path.splitext(os.path.basename(gff_path))[0]
        fasta_path = os.path.splitext(gff_path)[0] + ".fasta"
        if os.path.isfile(fasta_path):
            ref_index[ref_name] = (fasta_path, gff_path)
        else:
            print(f"[aa_report] AVERTISMENT: nu găsesc {fasta_path} pentru {gff_path}, ignorat.")

    # Probe
    samples = find_samples(folder)
    if not samples:
        print(f"[aa_report] EROARE: nicio probă în {folder}/snippy/")
        raise SystemExit(1)
    print(f"[aa_report] Probe: {samples}")

    # Citește best_ref pentru fiecare probă
    sample_best_refs = read_sample_best_refs(folder, samples)
    print(f"[aa_report] Best refs per probă: {sample_best_refs}")

    # Grupează probele după referința lor
    from collections import defaultdict
    groups = defaultdict(list)
    for s in samples:
        groups[sample_best_refs.get(s)].append(s)

    # Citește genomurile per probă (o singură dată)
    sample_snippy = {}
    sample_classic = {}
    for s in samples:
        sample_snippy[s] = parse_fasta(
            os.path.join(folder, "snippy", s, "snps.consensus.fa")
        )
        classic_path = os.path.join(folder, "assembly", f"{s}.fasta")
        sample_classic[s] = parse_fasta(classic_path) if os.path.isfile(classic_path) else {}

    # Generează câte o secțiune per grup referință
    all_sections = []
    hiv_resistance_page = None      # va fi setat dacă se detectează HIV1

    for ref_name, group_samples in sorted(groups.items(), key=lambda x: (x[0] is None, x[0])):
        if ref_name is None:
            print(f"[aa_report] Probe fără best_ref: {group_samples} — ignorate.")
            continue
        if ref_name not in ref_index:
            print(f"[aa_report] Referința '{ref_name}' nu are .gff — "
                  f"probele {group_samples} nu vor apărea în aa_report.")
            continue

        fasta_path, gff_path = ref_index[ref_name]
        ref_genome = parse_fasta(fasta_path)
        features   = parse_gff_cds(gff_path)
        ref_display = os.path.basename(fasta_path)
        print(f"[aa_report] Referință: {fasta_path} | Probe: {group_samples}")
        print(f"[aa_report] Gene CDS: {[f['name'] for f in features]}")

        # Secțiunea de aminoacizi (fără Stanford)
        section_body = build_report(
            group_samples, ref_genome, sample_snippy, sample_classic,
            features, ref_display, full_doc=False
        )
        all_sections.append(section_body)

        # Detecție patogen — separat per grup
        pathogen_id, pathogen_display, pathogen_reason, pathogen_conf = detect_pathogen(
            ref_display, features, ref_genome
        )
        if args.force_hiv and pathogen_id != "HIV1":
            pathogen_id      = "HIV1"
            pathogen_conf    = "low"
            pathogen_display = PATHOGEN_DB["HIV1"]["display_name"]
            pathogen_reason  = "[DEMO] --force-hiv activat manual. " + (pathogen_reason or "")
            print("[aa_report] DEMO: --force-hiv activ, se forțează analiza HIV1.")
        print(f"[aa_report] Patogen detectat: {pathogen_id} ({pathogen_conf})")
        print(f"[aa_report] Motiv: {pathogen_reason}")

        # Generează pagina separată de rezistență dacă e HIV1
        if pathogen_id == "HIV1":
            resistance_results = check_resistance(
                pathogen_id, group_samples, sample_snippy, sample_classic, features
            )
            hiv_resistance_page = build_resistance_page(
                pathogen_id, resistance_results, pathogen_reason, pathogen_conf
            )

    if not all_sections:
        print("[aa_report] EROARE: nicio secțiune generată (niciun .gff potrivit pentru probele existente).")
        raise SystemExit(1)

    # Adaugă link către pagina HIV în aa_report.html dacă există
    hiv_link = ""
    if hiv_resistance_page:
        hiv_link = (
            '<div style="background:#eaf4fb;border-left:4px solid #3498db;'
            'padding:10px 16px;border-radius:4px;margin:20px 0;font-size:.95em">'
            '&#128302; <b>Analiză rezistență HIV-1 (Stanford HIVdb):</b> '
            '<a href="hiv_resistance_report.html" style="color:#2980b9;font-weight:bold">'
            'Deschide raportul de rezistență HIV &rarr;</a>'
            '</div>'
        )

    separator = '<hr style="margin-top:60px;border-top:3px solid #2c3e50">'
    combined_body = separator.join(all_sections)

    full_html = f"""<!DOCTYPE html>
<html lang="ro">
<head>
<meta charset="UTF-8">
<title>VeGAS — Aminoacizi REF vs Snippy vs Clasic</title>
{CSS}
</head>
<body>
{hiv_link}
{combined_body}
</body>
</html>
"""

    # Scrie aa_report.html
    out_path = os.path.join(folder, "aa_report.html")
    with open(out_path, "w") as fh:
        fh.write(full_html)
    print(f"[aa_report] Scris: {out_path}")

    # Scrie hiv_resistance_report.html dacă e cazul
    if hiv_resistance_page:
        hiv_out = os.path.join(folder, "hiv_resistance_report.html")
        with open(hiv_out, "w") as fh:
            fh.write(hiv_resistance_page)
        print(f"[aa_report] Scris: {hiv_out}")


if __name__ == "__main__":
    main()
