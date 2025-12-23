#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
from dataclasses import dataclass
from typing import Dict, List

# -------------------------
# FASTA utilities (no deps)
# -------------------------
def read_fasta_lengths(path: str) -> Dict[str, int]:
    """Return {seqid: length} using first token after '>' as seqid."""
    lengths: Dict[str, int] = {}
    if not path or (not os.path.isfile(path)) or os.path.getsize(path) == 0:
        return lengths

    seqid = None
    seqlen = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seqid is not None:
                    lengths[seqid] = seqlen
                header = line[1:].strip()
                seqid = header.split()[0]
                seqlen = 0
            else:
                seqlen += len(line)
        if seqid is not None:
            lengths[seqid] = seqlen
    return lengths

# -------------------------
# Parse domtblout (HMMER)
# -------------------------
@dataclass
class DomHit:
    target_name: str
    target_acc: str
    query_name: str
    query_acc: str
    full_evalue: float
    full_score: float
    full_bias: float
    dom_ievalue: float
    dom_score: float
    dom_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    acc: float
    desc: str

def parse_domtblout(path: str) -> List[DomHit]:
    hits: List[DomHit] = []
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return hits

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 22:
                continue

            try:
                target_name = parts[0]
                target_acc = parts[1]
                query_name = parts[3]
                query_acc = parts[4]
                full_evalue = float(parts[6].replace("E", "e"))
                full_score = float(parts[7])
                full_bias = float(parts[8])
                dom_ievalue = float(parts[12].replace("E", "e"))
                dom_score = float(parts[13])
                dom_bias = float(parts[14])
                hmm_from = int(parts[15]); hmm_to = int(parts[16])
                ali_from = int(parts[17]); ali_to = int(parts[18])
                env_from = int(parts[19]); env_to = int(parts[20])
                acc = float(parts[21])
                desc = " ".join(parts[22:]) if len(parts) > 22 else ""
            except Exception:
                continue

            hits.append(DomHit(
                target_name=target_name, target_acc=target_acc,
                query_name=query_name, query_acc=query_acc,
                full_evalue=full_evalue, full_score=full_score, full_bias=full_bias,
                dom_ievalue=dom_ievalue, dom_score=dom_score, dom_bias=dom_bias,
                hmm_from=hmm_from, hmm_to=hmm_to, ali_from=ali_from, ali_to=ali_to,
                env_from=env_from, env_to=env_to, acc=acc, desc=desc
            ))
    return hits

def best_hmmsearch_per_seq(domtbl_path: str) -> Dict[str, DomHit]:
    """Return best domain hit per target sequence (min dom_ievalue, tie -> max dom_score)."""
    best: Dict[str, DomHit] = {}
    for h in parse_domtblout(domtbl_path):
        key = h.target_name
        if key not in best:
            best[key] = h
        else:
            cur = best[key]
            if (h.dom_ievalue < cur.dom_ievalue) or (h.dom_ievalue == cur.dom_ievalue and h.dom_score > cur.dom_score):
                best[key] = h
    return best

def best_pfam_for_accession(domtbl_path: str, pfam_id: str) -> Dict[str, DomHit]:
    """
    Return best Pfam hit for a given PFAM accession per query sequence.
    In Pfam hmmscan domtblout: target = Pfam family, query = your protein.
    We'll match target accession like PF01083.23 -> PF01083
    """
    pfam_id = pfam_id.strip()
    best: Dict[str, DomHit] = {}
    for h in parse_domtblout(domtbl_path):
        acc = re.sub(r"\..*$", "", h.target_acc)  # PF01083.23 -> PF01083
        if acc != pfam_id:
            continue
        q = h.query_name  # query is protein id
        if q not in best:
            best[q] = h
        else:
            cur = best[q]
            if (h.dom_ievalue < cur.dom_ievalue) or (h.dom_ievalue == cur.dom_ievalue and h.dom_score > cur.dom_score):
                best[q] = h
    return best

# -------------------------
# SignalP summary parser
# -------------------------
@dataclass
class SignalPSummary:
    sp_start: int
    sp_end: int
    cleavage_after_aa: int
    original_len: int
    new_len: int

def parse_signalp_summary(tsv_path: str) -> Dict[str, SignalPSummary]:
    m: Dict[str, SignalPSummary] = {}
    if not os.path.isfile(tsv_path) or os.path.getsize(tsv_path) == 0:
        return m
    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        _header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue

            # IMPORTANT: ensure seqid matches FASTA first-token IDs
            sid = parts[0].split()[0]

            try:
                sp_start = int(parts[1])
                sp_end = int(parts[2])
                cleave = int(parts[3])
                olen = int(parts[4])
                nlen = int(parts[5])
            except Exception:
                continue

            m[sid] = SignalPSummary(sp_start, sp_end, cleave, olen, nlen)
    return m

# -------------------------
# Main
# -------------------------
def find_repo_root() -> str:
    here = os.getcwd()
    for _ in range(6):
        if os.path.isdir(os.path.join(here, "data")) and os.path.isdir(os.path.join(here, "results")):
            return here
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return os.getcwd()

def main():
    ap = argparse.ArgumentParser(description="Build a CSV summary table for selected cutinase candidates.")
    ap.add_argument("--pfam", required=True, help="Pfam accession to select (e.g., PF01083)")
    ap.add_argument("--out", default="results/summary/cutinase_candidates_summary.csv", help="Output CSV path")
    args = ap.parse_args()

    root = find_repo_root()
    out_csv = os.path.join(root, args.out)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    # Samples determined by existence of final selected FASTA
    selected_glob = os.path.join(root, "results", "03_pfam", "*", "*_pfam_filtered.fasta")
    selected_fastas = sorted(glob.glob(selected_glob))

    if not selected_fastas:
        raise SystemExit(
            f"ERROR: No encuentro FASTAs finales en {selected_glob}\n"
            f"¿Ya corriste etapa 03 (Pfam) y se generó <sample>_pfam_filtered.fasta?"
        )

    rows: List[Dict[str, object]] = []

    for sel_fa in selected_fastas:
        sample = os.path.basename(os.path.dirname(sel_fa))

        # Per-sample paths
        hmm_domtbl = os.path.join(root, "results", "01_hmmsearch", sample, f"{sample}.domtblout")
        signalp_tsv = os.path.join(root, "results", "02_signalp", sample, f"{sample}_signalP_summary.tsv")
        trimmed_fa = os.path.join(root, "results", "02_signalp", sample, f"{sample}_signalp_trimmed.fasta")
        pfam_domtbl = os.path.join(root, "results", "03_pfam", sample, f"{sample}_pfam.domtblout")

        # Original proteome (optional)
        orig_prot = None
        for ext in ("faa", "fasta", "fa", "fna"):
            cand = os.path.join(root, "data", "proteomes", f"{sample}.{ext}")
            if os.path.isfile(cand):
                orig_prot = cand
                break

        # Load maps
        selected_len = read_fasta_lengths(sel_fa)
        selected_ids = list(selected_len.keys())

        trimmed_len = read_fasta_lengths(trimmed_fa)
        original_len = read_fasta_lengths(orig_prot) if orig_prot else {}

        sig = parse_signalp_summary(signalp_tsv)
        hmm_best = best_hmmsearch_per_seq(hmm_domtbl) if os.path.isfile(hmm_domtbl) else {}
        pf_best = best_pfam_for_accession(pfam_domtbl, args.pfam) if os.path.isfile(pfam_domtbl) else {}

        for sid in selected_ids:
            hmm = hmm_best.get(sid)
            pf = pf_best.get(sid)
            sp = sig.get(sid)

            # Fill original length: prefer original proteome, else use SignalP TSV
            orig_len_val = original_len.get(sid, "")
            if (orig_len_val == "" or orig_len_val is None) and sp is not None:
                orig_len_val = sp.original_len

            row = {
                "sample": sample,
                "sequence_id": sid,

                # lengths
                "length_original_aa": orig_len_val,
                "length_secreted_trimmed_aa": trimmed_len.get(sid, ""),
                "length_final_selected_aa": selected_len.get(sid, ""),

                # HMMER (hmmsearch)
                "hmm_query": getattr(hmm, "query_name", ""),
                "hmm_full_evalue": getattr(hmm, "full_evalue", ""),
                "hmm_full_bitscore": getattr(hmm, "full_score", ""),
                "hmm_domain_ievalue": getattr(hmm, "dom_ievalue", ""),
                "hmm_domain_bitscore": getattr(hmm, "dom_score", ""),
                "hmm_ali_from": getattr(hmm, "ali_from", ""),
                "hmm_ali_to": getattr(hmm, "ali_to", ""),
                "hmm_hmm_from": getattr(hmm, "hmm_from", ""),
                "hmm_hmm_to": getattr(hmm, "hmm_to", ""),
                "hmm_acc": getattr(hmm, "acc", ""),

                # SignalP
                "has_signal_peptide": "yes" if sp is not None else "no",
                "signalp_start": getattr(sp, "sp_start", ""),
                "signalp_end": getattr(sp, "sp_end", ""),
                "signalp_cleavage_after_aa": getattr(sp, "cleavage_after_aa", ""),
                "signalp_mature_length_aa": getattr(sp, "new_len", ""),

                # Pfam (cutinase domain)
                "pfam_accession": args.pfam,
                "pfam_hit": getattr(pf, "target_name", ""),
                "pfam_evalue_full": getattr(pf, "full_evalue", ""),
                "pfam_bitscore_full": getattr(pf, "full_score", ""),
                "pfam_domain_ievalue": getattr(pf, "dom_ievalue", ""),
                "pfam_domain_bitscore": getattr(pf, "dom_score", ""),
                "pfam_ali_from": getattr(pf, "ali_from", ""),
                "pfam_ali_to": getattr(pf, "ali_to", ""),
                "pfam_desc": getattr(pf, "desc", ""),
            }
            rows.append(row)

    fieldnames = [
        "sample","sequence_id",
        "length_original_aa","length_secreted_trimmed_aa","length_final_selected_aa",

        "hmm_query","hmm_full_evalue","hmm_full_bitscore","hmm_domain_ievalue","hmm_domain_bitscore",
        "hmm_ali_from","hmm_ali_to","hmm_hmm_from","hmm_hmm_to","hmm_acc",

        "has_signal_peptide","signalp_start","signalp_end","signalp_cleavage_after_aa","signalp_mature_length_aa",

        "pfam_accession","pfam_hit","pfam_evalue_full","pfam_bitscore_full","pfam_domain_ievalue","pfam_domain_bitscore",
        "pfam_ali_from","pfam_ali_to","pfam_desc",
    ]

    with open(out_csv, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"OK: wrote {len(rows)} rows -> {out_csv}")

if __name__ == "__main__":
    main()
