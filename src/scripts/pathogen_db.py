#!/usr/bin/env python3
"""
pathogen_db.py — Baza de date a mutațiilor de rezistență pentru multiple patogeni.

Structura PATHOGEN_DB[pathogen_id]:
  display_name     : Numele complet afișat în raport
  short_name       : Abrevierea afișată în badge-ul de confirmare
  keywords         : Liste de cuvinte-cheie (UPPERCASE) pentru detectare automată
                     în numele fișierului de referință sau în header-ul FASTA
  genome_size      : (min_bp, max_bp) sau None dacă nu e aplicabil
  gene_aliases     : {alias_UPPERCASE → canonical_name}  — mapare gene GFF → canonical
  required_genes   : Gene canonice necesare pentru confirmare
  drug_classes     : {cod → {"bg":str, "fg":str, "name":str}}
  drug_class_order : Ordinea de afișare a claselor de droguri
  source           : Sursa bazei de date (text afișat în raport)
  source_url       : URL sursă (link afișat în raport)
  mutations        : [(gene, pos_1based, wt_AA, alts_slash_sep,
                       drug_class_code, [drugs_list], level, note), ...]

Nivel: "major"    = rezistență clinică semnificativă
       "moderate" = rezistență intermediară / depend de context
       "minor"    = contribuție parțială, de obicei în combinație

Sursă: Stanford HIVdb (hivdb.stanford.edu) — API GraphQL confirmat.
  Date descărcate automat via fetch_resistance_db.py.
"""

PATHOGEN_DB = {}

# ═══════════════════════════════════════════════════════════════════════════════
# HIV-1
# ═══════════════════════════════════════════════════════════════════════════════
PATHOGEN_DB["HIV1"] = {
    "display_name": "HIV-1 (Human Immunodeficiency Virus type 1)",
    "short_name":   "HIV-1",
    "keywords": [
        "HIV", "HIV1", "HIV-1", "HIV_1", "HXB2", "HXB-2",
        "K03455", "AF033819", "U26942",
        "NL43", "NL4-3", "LAI", "IIIB", "BH10",
        "HUMAN IMMUNODEFICIENCY",
    ],
    "genome_size": (8500, 10500),
    "gene_aliases": {
        "PR": "PR", "PRO": "PR", "PROTEASE": "PR",
        "RT": "RT", "REVERSE_TRANSCRIPTASE": "RT", "REVERSE TRANSCRIPTASE": "RT",
        "IN": "IN", "INT": "IN", "INTEGRASE": "IN",
    },
    "required_genes": ["PR", "RT", "IN"],
    "drug_classes": {
        "PI":    {"bg": "#c0392b", "fg": "#fff", "name": "Inhibitori de Protează (PI)"},
        "NRTI":  {"bg": "#1a5276", "fg": "#fff", "name": "Inhibitori Nucleozidici RT (NRTI)"},
        "NNRTI": {"bg": "#6c3483", "fg": "#fff", "name": "Inhibitori Non-Nucleozidici RT (NNRTI)"},
        "INSTI": {"bg": "#1e8449", "fg": "#fff", "name": "Inhibitori ai Integrazei (INSTI)"},
    },
    "drug_class_order": ["PI", "NRTI", "NNRTI", "INSTI"],
    "source":     "Stanford HIVdb HIVDB_10.2 (publicat 2026-04-26)",
    "source_url": "https://hivdb.stanford.edu",
    "mutations": [
        # ── PR ──────────────────────────────────────────────
        ("PR",  30, "D", "N", "PI", ['NFV'], "major", ""),
        ("PR",  32, "V", "I", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV'], "moderate", ""),
        ("PR",  33, "L", "F", "PI", ['FPV/r', 'NFV', 'TPV/r'], "minor", ""),
        ("PR",  46, "M", "I", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r'], "moderate", ""),
        ("PR",  46, "M", "L", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "minor", ""),
        ("PR",  47, "I", "A", "PI", ['DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'TPV/r'], "major", ""),
        ("PR",  47, "I", "V", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'TPV/r'], "moderate", ""),
        ("PR",  48, "G", "V", "PI", ['ATV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r'], "major", ""),
        ("PR",  50, "I", "L", "PI", ['ATV/r'], "major", ""),
        ("PR",  50, "I", "V", "PI", ['DRV/r', 'FPV/r', 'LPV/r', 'NFV', 'SQV/r'], "major", ""),
        ("PR",  54, "I", "L", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r'], "major", ""),
        ("PR",  54, "I", "M", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "major", ""),
        ("PR",  58, "Q", "E", "PI", ['NFV', 'TPV/r'], "minor", ""),
        ("PR",  74, "T", "P", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'NFV', 'SQV/r', 'TPV/r'], "minor", ""),
        ("PR",  76, "L", "V", "PI", ['DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV'], "major", ""),
        ("PR",  82, "V", "A", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r'], "moderate", ""),
        ("PR",  82, "V", "F", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r'], "moderate", ""),
        ("PR",  82, "V", "L", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "moderate", ""),
        ("PR",  82, "V", "T", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "moderate", ""),
        ("PR",  82, "V", "S", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "moderate", ""),
        ("PR",  83, "N", "D", "PI", ['ATV/r', 'IDV/r', 'NFV', 'SQV/r', 'TPV/r'], "minor", ""),
        ("PR",  84, "I", "V", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "major", ""),
        ("PR",  84, "I", "A", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "major", ""),
        ("PR",  84, "I", "C", "PI", ['ATV/r', 'DRV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r', 'TPV/r'], "major", ""),
        ("PR",  88, "N", "D", "PI", ['ATV/r', 'NFV', 'SQV/r'], "major", ""),
        ("PR",  88, "N", "S", "PI", ['ATV/r', 'IDV/r', 'NFV', 'SQV/r'], "major", ""),
        ("PR",  90, "L", "M", "PI", ['ATV/r', 'FPV/r', 'IDV/r', 'LPV/r', 'NFV', 'SQV/r'], "major", ""),
        # ── RT ──────────────────────────────────────────────
        ("RT",  41, "M", "L", "NRTI", ['AZT', 'D4T', 'DDI'], "minor", ""),
        ("RT",  65, "K", "R", "NRTI", ['ABC', 'D4T', 'DDI', 'FTC', '3TC', 'TDF'], "major", ""),
        ("RT",  65, "K", "E", "NRTI", ['ABC', 'D4T', 'DDI', 'TDF'], "minor", ""),
        ("RT",  65, "K", "N", "NRTI", ['ABC', 'D4T', 'DDI', 'FTC', '3TC', 'TDF'], "moderate", ""),
        ("RT",  67, "D", "N", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT",  67, "D", "E", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT",  67, "D", "G", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT",  67, "D", "S", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT",  70, "K", "R", "NRTI", ['AZT', 'D4T', 'DDI'], "moderate", ""),
        ("RT",  70, "K", "E", "NRTI", ['ABC', 'D4T', 'DDI', 'FTC', '3TC', 'TDF'], "minor", ""),
        ("RT",  74, "L", "V", "NRTI", ['ABC', 'DDI'], "major", ""),
        ("RT",  74, "L", "I", "NRTI", ['ABC', 'DDI'], "major", ""),
        ("RT",  77, "F", "L", "NRTI", ['D4T', 'DDI'], "minor", ""),
        ("RT", 115, "Y", "F", "NRTI", ['ABC', 'TDF'], "moderate", ""),
        ("RT", 116, "F", "Y", "NRTI", ['AZT', 'D4T', 'DDI'], "minor", ""),
        ("RT", 151, "Q", "M", "NRTI", ['ABC', 'AZT', 'D4T', 'DDI', 'FTC', '3TC', 'TDF'], "major", ""),
        ("RT", 184, "M", "V", "NRTI", ['ABC', 'DDI', 'FTC', 'ISL', '3TC'], "major", ""),
        ("RT", 184, "M", "I", "NRTI", ['ABC', 'DDI', 'FTC', 'ISL', '3TC'], "major", ""),
        ("RT", 210, "L", "W", "NRTI", ['AZT', 'D4T', 'DDI'], "minor", ""),
        ("RT", 215, "T", "Y", "NRTI", ['ABC', 'AZT', 'D4T', 'DDI', 'ISL', 'TDF'], "major", ""),
        ("RT", 215, "T", "F", "NRTI", ['ABC', 'AZT', 'D4T', 'DDI', 'ISL', 'TDF'], "major", ""),
        ("RT", 219, "K", "Q", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT", 219, "K", "E", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT", 219, "K", "N", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT", 219, "K", "R", "NRTI", ['AZT', 'D4T'], "minor", ""),
        ("RT", 100, "K", "I", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 101, "K", "E", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "moderate", ""),
        ("RT", 101, "K", "P", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 103, "K", "N", "NNRTI", ['EFV', 'NVP'], "major", ""),
        ("RT", 103, "K", "S", "NNRTI", ['EFV', 'NVP'], "major", ""),
        ("RT", 106, "V", "A", "NNRTI", ['DOR', 'EFV', 'NVP'], "major", ""),
        ("RT", 106, "V", "M", "NNRTI", ['DOR', 'DPV', 'EFV', 'NVP'], "major", ""),
        ("RT", 108, "V", "I", "NNRTI", ['EFV', 'NVP'], "minor", ""),
        ("RT", 138, "E", "A", "NNRTI", ['DPV', 'ETR', 'RPV'], "minor", ""),
        ("RT", 138, "E", "G", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 138, "E", "K", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "moderate", ""),
        ("RT", 138, "E", "Q", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 138, "E", "R", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 179, "V", "D", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 179, "V", "F", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 179, "V", "L", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 181, "Y", "C", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 181, "Y", "I", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 181, "Y", "V", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 188, "Y", "C", "NNRTI", ['EFV', 'NVP'], "major", ""),
        ("RT", 188, "Y", "L", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 188, "Y", "H", "NNRTI", ['EFV', 'NVP'], "major", ""),
        ("RT", 190, "G", "A", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 190, "G", "E", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 190, "G", "S", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 190, "G", "Q", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 221, "H", "Y", "NNRTI", ['DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "minor", ""),
        ("RT", 225, "P", "H", "NNRTI", ['DOR', 'DPV', 'EFV', 'NVP'], "moderate", ""),
        ("RT", 227, "F", "C", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 227, "F", "L", "NNRTI", ['DOR', 'DPV', 'EFV', 'NVP'], "major", ""),
        ("RT", 230, "M", "I", "NNRTI", ['DPV', 'EFV', 'NVP', 'RPV'], "moderate", ""),
        ("RT", 230, "M", "L", "NNRTI", ['DOR', 'DPV', 'EFV', 'ETR', 'NVP', 'RPV'], "major", ""),
        ("RT", 238, "K", "T", "NNRTI", ['EFV', 'NVP'], "moderate", ""),
        ("RT", 318, "Y", "F", "NNRTI", ['DOR', 'EFV', 'NVP'], "major", ""),
        # ── IN ──────────────────────────────────────────────
        ("IN",  66, "T", "I", "INSTI", ['CAB', 'EVG', 'RAL'], "major", ""),
        ("IN",  66, "T", "A", "INSTI", ['EVG', 'RAL'], "major", ""),
        ("IN",  66, "T", "K", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN",  92, "E", "Q", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN",  92, "E", "G", "INSTI", ['CAB', 'EVG', 'RAL'], "moderate", ""),
        ("IN", 121, "F", "Y", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN", 140, "G", "S", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "moderate", ""),
        ("IN", 140, "G", "A", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "moderate", ""),
        ("IN", 140, "G", "C", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "moderate", ""),
        ("IN", 143, "Y", "R", "INSTI", ['EVG', 'RAL'], "major", ""),
        ("IN", 143, "Y", "C", "INSTI", ['EVG', 'RAL'], "major", ""),
        ("IN", 143, "Y", "H", "INSTI", ['EVG', 'RAL'], "major", ""),
        ("IN", 147, "S", "G", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN", 148, "Q", "H", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN", 148, "Q", "R", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN", 148, "Q", "K", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN", 155, "N", "H", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
        ("IN", 155, "N", "S", "INSTI", ['CAB', 'EVG', 'RAL'], "moderate", ""),
        ("IN", 157, "E", "Q", "INSTI", ['EVG', 'RAL'], "minor", ""),
        ("IN", 263, "R", "K", "INSTI", ['BIC', 'CAB', 'DTG', 'EVG', 'RAL'], "major", ""),
    ],
}

# ═══════════════════════════════════════════════════════════════════════════════
# Detecție automată
# ═══════════════════════════════════════════════════════════════════════════════

def detect_pathogen(ref_name, features, ref_genome):
    """Detectează automat patogenul pe baza referinței, genelor GFF și genomului.

    Returnează: (pathogen_id, display_name, reason, confidence)
      pathogen_id  : cheie din PATHOGEN_DB sau None
      display_name : numele complet sau None
      reason       : string cu semnalele găsite
      confidence   : 'high', 'low' sau 'none'
    """
    upper_ref     = ref_name.upper()
    fasta_headers = " ".join(ref_genome.keys()).upper()
    gene_names    = {f["name"].upper() for f in features}

    scores = {}   # pathogen_id → (punctaj, [semnale])

    for pid, pdata in PATHOGEN_DB.items():
        signals = []

        # 1. Keyword în numele fișierului referință
        for kw in pdata["keywords"]:
            if kw.upper() in upper_ref:
                signals.append(f"Keyword '{kw}' în numele referinței")
                break

        # 2. Keyword în header-ul FASTA
        for kw in pdata["keywords"]:
            if kw.upper() in fasta_headers:
                signals.append(f"Keyword '{kw}' în header FASTA")
                break

        # 3. Gene recunoscute în GFF
        aliases = pdata["gene_aliases"]
        found   = list({aliases[g] for g in gene_names if g in aliases})
        if found:
            signals.append(f"Gene {found} în GFF")

        # 4. Dimensiunea genomului
        if pdata["genome_size"]:
            gmin, gmax = pdata["genome_size"]
            for seq in ref_genome.values():
                if gmin <= len(seq) <= gmax:
                    signals.append(
                        f"Lungimea genomului {len(seq)} bp în intervalul {gmin}–{gmax} bp"
                    )
                    break

        if signals:
            # Semnal keyword = prioritate mare; fiecare semnal = +1 punct
            kw_match = any("Keyword" in s for s in signals)
            score    = len(signals) + (3 if kw_match else 0)
            scores[pid] = (score, signals)

    if not scores:
        return None, None, "Niciun patogen cunoscut detectat în referință.", "none"

    best          = max(scores, key=lambda p: scores[p][0])
    score, signals = scores[best]
    reason        = "; ".join(signals)
    display       = PATHOGEN_DB[best]["display_name"]
    confidence    = "high" if score >= 4 else "low"

    return best, display, reason, confidence
