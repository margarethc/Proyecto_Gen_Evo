#!/usr/bin/env python3
"""
extract_from_list.py

Extract sequences from a FASTA (.faa) file given a list of sequence names
and write them to an output FASTA.

Features:
- Provide names either as a plain text file (--list-file) with one name per line,
  or directly via the command line (--names "id1,id2,id3").
- Matching is by default exact against the first token of the FASTA header
  (the part before the first whitespace). Options available for substring or
  case-insensitive matching.
- Output preserves the order of names supplied (duplicates in the input list
  are ignored, but order is kept).
- If you prefer to preserve FASTA input order rather than list order, use --preserve-fasta-order.

Usage examples:
    # Using a file with one id per line
    python3 extract_from_list.py --list-file ids.txt --faa sequences.faa -o extracted.faa

    # Using names on the command line
    python3 extract_from_list.py --names "LMB_CAL9_8192.t1,LMB_CAL9_7202.t1" --faa sequences.faa

    # Case-insensitive substring matching, output in FASTA order
    python3 extract_from_list.py --list-file ids.txt --faa sequences.faa -o out.faa --match-substring --ignore-case --preserve-fasta-order
"""
from argparse import ArgumentParser
import sys

def read_name_list_from_file(path):
    names = []
    with open(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            names.append(line)
    return names

def parse_names_arg(names_arg):
    # comma or whitespace separated
    parts = []
    for chunk in names_arg.split(','):
        for token in chunk.split():
            t = token.strip()
            if t:
                parts.append(t)
    return parts

def fasta_iter(path):
    """Yield (header_without_gt, sequence) for each FASTA record. Sequence is a single-line string."""
    with open(path, 'rt') as fh:
        header = None
        seq_lines = []
        for raw in fh:
            if raw.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = raw[1:].rstrip('\n')
                seq_lines = []
            else:
                seq_lines.append(raw.strip())
        if header is not None:
            yield header, ''.join(seq_lines)

def index_fasta_by_id(faa_path, id_selector=lambda header: header.split()[0]):
    """
    Build a dict mapping id_selector(header) -> (header, seq)
    Note: this will load matched sequences into memory. For very large FASTA files
    consider streaming methods.
    """
    idx = {}
    for header, seq in fasta_iter(faa_path):
        key = id_selector(header)
        if key not in idx:
            idx[key] = (header, seq)
    return idx

def write_fasta_record(outfh, header, seq):
    outfh.write(f">{header}\n")
    outfh.write(seq + "\n")

def main():
    p = ArgumentParser(description="Extract sequences from a .faa given a list of names")
    p.add_argument('--list-file', help="Plain text file with one sequence name per line (optional)")
    p.add_argument('--names', help="Comma- or space-separated list of names (optional), e.g. 'id1,id2'")
    p.add_argument('--faa', required=True, help="Input FASTA (.faa) file")
    p.add_argument('-o', '--out', default='extracted.faa', help="Output FASTA file")
    p.add_argument('--match-substring', action='store_true', help="Match names as substring anywhere in header (slower)")
    p.add_argument('--ignore-case', action='store_true', help="Case-insensitive matching")
    p.add_argument('--preserve-fasta-order', action='store_true',
                   help="Write matches in the order they appear in the FASTA (default is list order)")
    p.add_argument('--verbose', action='store_true', help="Print progress")
    args = p.parse_args()

    # Build target list
    targets = []
    if args.list_file:
        targets += read_name_list_from_file(args.list_file)
    if args.names:
        targets += parse_names_arg(args.names)
    if not targets:
        print("No target names provided. Use --list-file or --names.", file=sys.stderr)
        sys.exit(1)

    # Remove duplicates but preserve order
    seen = set()
    ordered_targets = []
    for t in targets:
        if t not in seen:
            ordered_targets.append(t)
            seen.add(t)

    if args.ignore_case:
        # Work in lowercase for matching
        ordered_targets_keyed = [t.lower() for t in ordered_targets]
    else:
        ordered_targets_keyed = ordered_targets

    if args.verbose:
        print(f"Targets: {len(ordered_targets)} unique names", file=sys.stderr)

    # If preserving FASTA order, stream through FASTA and write records as found.
    if args.preserve_fasta_order:
        matched = set()
        with open(args.out, 'wt') as outfh:
            for header, seq in fasta_iter(args.faa):
                hay = header.lower() if args.ignore_case else header
                header_id = header.split()[0]
                hay_id = header_id.lower() if args.ignore_case else header_id

                is_match = False
                if args.match_substring:
                    # any target substring in full header
                    for tkey in ordered_targets_keyed:
                        if tkey in hay:
                            is_match = True
                            matched.add(tkey)
                            break
                else:
                    # exact match against first token
                    if hay_id in ordered_targets_keyed:
                        is_match = True
                        matched.add(hay_id)

                if is_match:
                    write_fasta_record(outfh, header, seq)

        if args.verbose:
            print(f"Wrote {len(matched)} sequences to {args.out} (FASTA order)", file=sys.stderr)
        return

    # Otherwise, prefer list order: index FASTA by id or by header and then write in list order.
    if args.match_substring:
        # We will index by full header (or lowercase header) and check substring matches.
        # For efficiency, build a dict mapping header_key -> (header, seq)
        if args.verbose:
            print("Indexing FASTA by full header (may use more memory)...", file=sys.stderr)
        if args.ignore_case:
            id_selector = lambda header: header.lower()
        else:
            id_selector = lambda header: header
        fasta_index = index_fasta_by_id(args.faa, id_selector=id_selector)

        written = 0
        with open(args.out, 'wt') as outfh:
            for tkey in ordered_targets_keyed:
                # find first header containing tkey
                found = False
                for header_key, (header, seq) in fasta_index.items():
                    if tkey in header_key:
                        write_fasta_record(outfh, header, seq)
                        written += 1
                        found = True
                        break
                if args.verbose and not found:
                    print(f"Warning: '{tkey}' not found in FASTA (substring search).", file=sys.stderr)

        if args.verbose:
            print(f"Wrote {written} sequences to {args.out} (list order, substring matching)", file=sys.stderr)
        return

    # Default exact-id matching: index by first token of header
    if args.verbose:
        print("Indexing FASTA by first token of header (id)...", file=sys.stderr)
    if args.ignore_case:
        id_selector = lambda header: header.split()[0].lower()
    else:
        id_selector = lambda header: header.split()[0]
    fasta_index = index_fasta_by_id(args.faa, id_selector=id_selector)

    written = 0
    with open(args.out, 'wt') as outfh:
        for orig_name, tkey in zip(ordered_targets, ordered_targets_keyed):
            if tkey in fasta_index:
                header, seq = fasta_index[tkey]
                write_fasta_record(outfh, header, seq)
                written += 1
            else:
                if args.verbose:
                    print(f"Warning: '{orig_name}' not found in FASTA (exact id match).", file=sys.stderr)

    if args.verbose:
        print(f"Wrote {written} sequences to {args.out} (list order, exact id matching)", file=sys.stderr)

if __name__ == '__main__':
    main()
