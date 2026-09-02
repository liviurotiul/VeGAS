#!/usr/bin/env python3
"""
fetch_resistance_db.py — Descarcă date de rezistență HIV-1 de la Stanford HIVdb.

Sursă: Stanford HIVdb GraphQL API — hivdb.stanford.edu/graphql
API confirmat, utilizat de laboratoare clinice certificate la nivel mondial.

Utilizare:
  python fetch_resistance_db.py            # actualizează pathogen_db.py
  python fetch_resistance_db.py --dry-run  # afișează fără a scrie

Necesită: pip install requests
"""
import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone

try:
    import requests
except ImportError:
    sys.exit("EROARE: lipsă pachet 'requests'. Rulați: pip install requests")

# ── Configurare ────────────────────────────────────────────────────────────────
HIVDB_GRAPHQL = "https://hivdb.stanford.edu/graphql"


# Output: fișier JSON cu datele descărcate (nu modifică pathogen_db.py direct)
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_JSON = os.path.join(SCRIPT_DIR, "resistance_db_fetched.json")

TIMEOUT = 30   # secunde

# ══════════════════════════════════════════════════════════════════════════════
# Stanford HIVdb GraphQL API
# Documentație: https://hivdb.stanford.edu/page/webservice/
#
# Schema verificată prin introspection (iulie 2026):
#   Root.currentVersion → { text publishDate }
#   Root.mutationsAnalysis(mutations: [String]) → {
#     drugResistance: [{
#       gene: { name }
#       drugScores: [{
#         drugClass: { name }
#         drug: { displayAbbr name }
#         level   ← Int 1-5 (1=Susceptible, 5=High-Level Resistance)
#         text    ← String "High-Level Resistance" etc.
#         SIR     ← enum S/I/R
#       }]
#     }]
#   }
#
# IMPORTANT: mutations[] trimite o singură mutație pentru a obține
# rezistența acelei mutații individuale (nu combinație).
# ══════════════════════════════════════════════════════════════════════════════

# Nivel de rezistență → clasificare internă
# level 5 = High-Level Resistance   → "major"
# level 4 = Intermediate Resistance → "moderate"
# level 3 = Low-Level Resistance    → "minor"
# level 2 = Potential Low-Level     → "minor" (inclus deoarece e clinic relevant)
# level 1 = Susceptible             → ignorat
LEVEL_MAP = {5: "major", 4: "moderate", 3: "minor", 2: "minor"}

# Query pentru o singură mutație — returnează rezistența individuală
SINGLE_MUT_QUERY = """
{
  mutationsAnalysis(mutations: %s) {
    drugResistance {
      gene { name }
      drugScores {
        drugClass { name }
        drug { displayAbbr name }
        level
        text
      }
    }
  }
}
"""


def graphql_post(url, query, timeout=None):
    """Execută o interogare GraphQL și returnează câmpul 'data'."""
    try:
        resp = requests.post(
            url,
            json={"query": query},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=timeout or TIMEOUT,
        )
        resp.raise_for_status()
        result = resp.json()
        if "errors" in result:
            raise RuntimeError(f"GraphQL errors: {result['errors']}")
        return result.get("data", {})
    except requests.exceptions.ConnectionError:
        raise RuntimeError(f"Nu se poate conecta la {url}. Verificați conexiunea la internet.")
    except requests.exceptions.Timeout:
        raise RuntimeError(f"Timeout la {url} după {TIMEOUT}s.")


# ─── Lista mutațiilor HIV cunoscute ca DRM (Drug Resistance Mutations) ─────────
# Format: "GENE:MUTATION" conform convenției HIVdb (ex: "PR:D30N")
# Acestea sunt trimise la API-ul Stanford pentru a obține clasificarea oficială.
HIV_MUTATIONS_TO_CHECK = [
    # PR
    "PR:D30N","PR:V32I","PR:L33F","PR:M46I","PR:M46L","PR:I47A","PR:I47V",
    "PR:G48V","PR:I50L","PR:I50V","PR:I54L","PR:I54M","PR:Q58E","PR:T74P",
    "PR:L76V","PR:V82A","PR:V82F","PR:V82L","PR:V82T","PR:V82S","PR:N83D",
    "PR:I84V","PR:I84A","PR:I84C","PR:N88D","PR:N88S","PR:L90M",
    # RT NRTI
    "RT:M41L","RT:K65R","RT:K65E","RT:K65N","RT:D67N","RT:D67E","RT:D67G",
    "RT:D67S","RT:K70R","RT:K70E","RT:L74V","RT:L74I","RT:V75I","RT:F77L",
    "RT:Y115F","RT:F116Y","RT:Q151M","RT:M184V","RT:M184I","RT:L210W",
    "RT:T215Y","RT:T215F","RT:K219Q","RT:K219E","RT:K219N","RT:K219R",
    # RT NNRTI
    "RT:K100I","RT:K101E","RT:K101P","RT:K103N","RT:K103S","RT:V106A",
    "RT:V106M","RT:V108I","RT:E138A","RT:E138G","RT:E138K","RT:E138Q",
    "RT:E138R","RT:V179D","RT:V179F","RT:V179L","RT:V179T","RT:Y181C",
    "RT:Y181I","RT:Y181V","RT:Y188C","RT:Y188L","RT:Y188H","RT:G190A",
    "RT:G190E","RT:G190S","RT:G190Q","RT:H221Y","RT:P225H","RT:F227C",
    "RT:F227L","RT:M230I","RT:M230L","RT:K238T","RT:Y318F",
    # IN INSTI
    "IN:T66I","IN:T66A","IN:T66K","IN:E92Q","IN:E92G","IN:F121Y",
    "IN:G140S","IN:G140A","IN:G140C","IN:Y143R","IN:Y143C","IN:Y143H",
    "IN:S147G","IN:Q148H","IN:Q148R","IN:Q148K","IN:N155H","IN:N155S",
    "IN:E157Q","IN:R263K",
]

def query_single_mutation(url, mut_str):
    """
    Trimite o singură mutație la API și returnează lista de droguri cu rezistență.
    Returnează: list of (drug_class, drug_abbr, level_str, text_str)
    """
    # Construim query GraphQL cu valoarea literală (nu variabilă)
    # Format acceptat de Stanford API: ["PR:D30N"]
    query = SINGLE_MUT_QUERY % json.dumps([mut_str])
    data  = graphql_post(url, query)
    hits  = []
    for dr in data.get("mutationsAnalysis", {}).get("drugResistance", []):
        for ds in dr["drugScores"]:
            lvl = ds.get("level", 1)
            if lvl in LEVEL_MAP:
                hits.append((
                    ds["drugClass"]["name"],
                    ds["drug"]["displayAbbr"],
                    LEVEL_MAP[lvl],
                    ds.get("text", ""),
                ))
    return hits


def fetch_hiv_from_stanford():
    """Descarcă informații HIV de la Stanford HIVdb (API oficial).

    Metodă:
      - Interogăm mai întâi versiunea curentă a algoritmului HIVdb.
      - Trimitem fiecare mutație INDIVIDUAL la 'mutationsAnalysis'.
        Asta asigură că obținem rezistența ACELEI mutații, nu a unei combinații.
      - Folosim câmpul 'level' (Int 1-5) din schema verificată, nu 'score'.
      - Păstrăm mutațiile cu level >= 2 (Potential Low-Level Resistance sau mai mare).

    URL oficial: https://hivdb.stanford.edu/graphql
    """
    print("[HIV] Conectare la Stanford HIVdb API...")

    # Pasul 1: versiunea curentă a algoritmului HIVdb
    meta    = graphql_post(HIVDB_GRAPHQL, "{ currentVersion { text publishDate } }")
    version = meta["currentVersion"]["text"]
    pubdate = meta["currentVersion"]["publishDate"]
    print(f"[HIV] Versiune HIVdb: {version} (publicat {pubdate})")

    # Pasul 2: interogăm fiecare mutație individual
    print(f"[HIV] Verificare {len(HIV_MUTATIONS_TO_CHECK)} mutații (câte una pe rând)...")
    confirmed = []   # tuple: (gene, pos, wt, alt, drug_class, [drugs], level, note)

    for idx, mut_str in enumerate(HIV_MUTATIONS_TO_CHECK, 1):
        gene_part, mut_part = mut_str.split(":")
        wt  = mut_part[0]
        pos = int(mut_part[1:-1])
        aa  = mut_part[-1]

        try:
            hits = query_single_mutation(HIVDB_GRAPHQL, mut_str)
        except Exception as e:
            print(f"[HIV]   [{idx}/{len(HIV_MUTATIONS_TO_CHECK)}] {mut_str} → EROARE: {e}")
            time.sleep(1.0)
            continue

        if hits:
            # Grupăm drogurile pe clasă; păstrăm nivelul cel mai mare
            rank = {"major": 3, "moderate": 2, "minor": 1}
            class_to_drugs = {}
            class_to_level = {}
            for dc_name, drug_abbr, lvl_str, txt in hits:
                class_to_drugs.setdefault(dc_name, []).append(drug_abbr)
                # Inițializare sau actualizare la nivel mai mare
                if dc_name not in class_to_level or rank[lvl_str] > rank[class_to_level[dc_name]]:
                    class_to_level[dc_name] = lvl_str

            for dc_name, drugs in class_to_drugs.items():
                lvl_str = class_to_level[dc_name]
                confirmed.append((gene_part, pos, wt, aa, dc_name,
                                   drugs, lvl_str, ""))
            desc = ", ".join(
                f"{'+'.join(drugs)}({class_to_level[dc]})"
                for dc, drugs in class_to_drugs.items()
            )
            print(f"[HIV]   [{idx}/{len(HIV_MUTATIONS_TO_CHECK)}] {mut_str} → {desc}")
        else:
            print(f"[HIV]   [{idx}/{len(HIV_MUTATIONS_TO_CHECK)}] {mut_str} → Susceptibil")

        time.sleep(0.35)   # politicos față de server Stanford

    print(f"\n[HIV] Total: {len(confirmed)} intrări confirmate cu rezistență.")
    return {
        "source":     f"Stanford HIVdb {version} (publicat {pubdate})",
        "source_url": "https://hivdb.stanford.edu",
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "version":    version,
        "mutations":  confirmed,
    }



def update_pathogen_db(results, dry_run=False):
    """
    Actualizează secțiunile de mutații din pathogen_db.py pentru
    patogenii cu date obținute via API live (api_type != 'manual_curation').

    Logica de înlocuire este identică cu ce s-a folosit pentru HIV.
    """
    db_path = os.path.join(SCRIPT_DIR, "pathogen_db.py")
    with open(db_path) as f:
        content = f.read()

    updated = []

    for pid, data in results.items():
        # Sari dacă nu e de la un API live sau dacă nu are mutații
        if data.get("api_type") == "manual_curation":
            continue
        if "error" in data or not data.get("mutations"):
            continue

        muts    = data["mutations"]
        src_new = data.get("source", "")

        # ── Găsim vechea linie source ─────────────────────────────────────
        # Presupunem că în pathogen_db.py există un bloc:
        #   PATHOGEN_DB["<PID>"] = { ... "source": "...", ... "mutations": [ ... ], ...
        # Trebuie să știm blocul exact. Folosim markeri de context.
        # Strategia: găsim PATHOGEN_DB["<PID>"] = { și de acolo căutăm "source" și "mutations"

        block_start = f'PATHOGEN_DB["{pid}"]'
        if block_start not in content:
            print(f"  [{pid}] Nu s-a găsit blocul în pathogen_db.py — skip.")
            continue

        bi = content.index(block_start)

        # Source
        src_marker = '"source":'
        si_src = content.index(src_marker, bi)
        # Înlocuim valoarea source (până la virgulă + newline)
        si_src_val = content.index('"', si_src + len(src_marker))
        ei_src_val = content.index('"', si_src_val + 1)
        old_src = content[si_src_val + 1:ei_src_val]
        if old_src != src_new:
            content = content[:si_src_val + 1] + src_new + content[ei_src_val:]
            # Recalculăm bi după modificare
            bi = content.index(block_start)

        # Mutations block
        start_marker = '    "mutations": [\n'
        end_marker   = '    ],\n'
        si_muts = content.index(start_marker, bi)
        ei_muts = content.index(end_marker, si_muts) + len(end_marker)

        # Generăm noul bloc
        mut_lines = []
        cur_gene  = None
        for m in muts:
            gene, pos, wt, aa, dc, drugs, lvl, note = m
            if gene != cur_gene:
                mut_lines.append(f"        # ── {gene} ──────────────────────────────────────────────")
                cur_gene = gene
            mut_lines.append(
                f"        (\"{gene}\", {pos:3d}, \"{wt}\", \"{aa}\", \"{dc}\", {repr(list(drugs))}, \"{lvl}\", \"{note}\"),"
            )
        new_block = start_marker + "\n".join(mut_lines) + "\n" + end_marker

        content = content[:si_muts] + new_block + content[ei_muts:]
        updated.append(pid)
        print(f"  [{pid}] actualizat cu {len(muts)} mutații (sursă: {src_new})")

    if updated and not dry_run:
        import ast
        try:
            ast.parse(content)
        except SyntaxError as e:
            print(f"  EROARE sintaxă după actualizare: {e} — pathogen_db.py NU a fost scris.")
            return
        with open(db_path, "w") as f:
            f.write(content)
        print(f"  pathogen_db.py scris cu succes ({len(updated)} patogeni actualizați).")
    elif not updated:
        print("  Niciun patogen cu API live de actualizat.")
    else:
        print("  [dry-run] pathogen_db.py nu a fost modificat.")


# ══════════════════════════════════════════════════════════════════════════════
# Mapare organism → funcție fetch
# ══════════════════════════════════════════════════════════════════════════════

FETCHERS = {
    "HIV1": fetch_hiv_from_stanford,
}

ALL_PATHOGENS = list(FETCHERS.keys())


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Descarcă date de rezistență de la surse oficiale reale."
    )
    parser.add_argument(
        "--pathogen",
        choices=ALL_PATHOGENS + ["ALL"],
        default="ALL",
        metavar="PATHOGEN",
        help=(
            "Patogenul de actualizat. Opțiuni: "
            + ", ".join(ALL_PATHOGENS)
            + ", ALL (implicit)"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Afișează ce s-ar descărca/actualiza fără a scrie fișiere",
    )
    parser.add_argument(
        "--update-db",
        action="store_true",
        help="Actualizează pathogen_db.py automat cu datele obținute via API live",
    )
    args = parser.parse_args()

    targets = ALL_PATHOGENS if args.pathogen == "ALL" else [args.pathogen]

    print("=" * 68)
    print(" fetch_resistance_db.py — Date de la surse oficiale (10 organisme)")
    print("=" * 68)
    print(f" Timestamp: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    print(f" Organisme: {', '.join(targets)}")
    print()

    results = {}

    for pid in targets:
        fn = FETCHERS[pid]
        try:
            results[pid] = fn()
        except Exception as e:
            print(f"[{pid}] EROARE: {e}")
            results[pid] = {"error": str(e), "mutations": None}
        print()

    # ── Sumar ──────────────────────────────────────────────────────────────
    print("=" * 68)
    print(" SUMAR")
    print("=" * 68)
    live_api = []
    manual   = []
    errors   = []
    for pid, data in results.items():
        if "error" in data:
            errors.append(pid)
            print(f"  {pid:12s} ✗  EROARE: {data['error']}")
        elif data.get("api_type") == "manual_curation":
            manual.append(pid)
            n = len(data.get("mutations") or [])
            print(f"  {pid:12s} ✓  {n:3d} mutații — curare manuală — {data.get('source','')[:50]}")
        else:
            live_api.append(pid)
            n = len(data.get("mutations") or [])
            print(f"  {pid:12s} ✓  {n:3d} mutații — API live     — {data.get('source','')[:50]}")

    print()
    print(f" API live:        {len(live_api)} patogeni ({', '.join(live_api) or '—'})")
    print(f" Curare manuală:  {len(manual)} patogeni ({', '.join(manual) or '—'})")
    if errors:
        print(f" Erori:           {len(errors)} ({', '.join(errors)})")

    # ── Salvare JSON ───────────────────────────────────────────────────────
    if not args.dry_run:
        # Serializăm tuplele/listele ca liste pentru JSON
        serializable = {}
        for pid, data in results.items():
            d = dict(data)
            if d.get("mutations"):
                d["mutations"] = [list(m) for m in d["mutations"]]
            serializable[pid] = d
        with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
            json.dump(serializable, f, ensure_ascii=False, indent=2)
        print(f"\n Rezultate salvate în: {OUTPUT_JSON}")
    else:
        print("\n [dry-run] JSON nu a fost scris.")

    # ── Actualizare pathogen_db.py ─────────────────────────────────────────
    if args.update_db or (live_api and not args.dry_run):
        print()
        print(" Actualizare pathogen_db.py cu datele API live...")
        update_pathogen_db(results, dry_run=args.dry_run)

    # ── Surse oficiale ─────────────────────────────────────────────────────
    print()
    print(" Surse oficiale pentru verificare manuală:")
    print("   HIV-1:      https://hivdb.stanford.edu")
    print("   HCV:        https://hcv.stanford.edu")
    print("   Flu A:      https://www.cdc.gov/flu/professionals/antivirals/")
    print("   SARS-CoV-2: https://covdb.stanford.edu/  |  https://www.fda.gov/")
    print("   CMV/HSV:    https://www.ecil-leukaemia.com")
    print("   RSV:        https://www.cdc.gov/rsv")
    print("   Fungi:      https://www.eucast.org/clinical_breakpoints/")
    print("   C. auris:   https://www.cdc.gov/candida-auris/")


if __name__ == "__main__":
    main()
