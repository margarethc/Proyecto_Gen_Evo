
import argparse
import sys

def fasta_iter(path):
    """Yield (header_without_gt, sequence) for each FASTA record."""
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

def parse_range(r):
    parts = r.split('-', 1)
    if len(parts) != 2:
        raise ValueError(f"Bad range format: '{r}' (expected A-B)")
    a = int(parts[0])
    b = int(parts[1])
    return a, b

def parse_positions_arg(pos_arg):
    """Parse a positions string like '32-33' or '10-20,40-50' into list of (a,b) tuples or None markers."""
    pos_arg = pos_arg.strip()
    if not pos_arg:
        return []
    ranges = []
    for piece in pos_arg.split(','):
        piece = piece.strip()
        if not piece:
            continue
        if piece == '-':
            ranges.append(None)
            continue
        ranges.append(parse_range(piece))
    return ranges

def read_positions_file(path):
    """
    Read a positions file. Recognizes three forms:
      - Per-id: "seq_id:32-33" or "seq_id:10-20,40-50"  -> per_seq dict
      - Ordered list (no ':' anywhere): each line applies to successive FASTA records.
        A line '-' becomes [None] meaning "no cut".
      - Global single line: if file has a single (non-empty, non-comment) line without ':'
        it will be treated as the global ranges.
    Returns a tuple: (mode, data)
      - mode == 'per_id' -> data is per_seq dict (seq_id -> ranges list)
      - mode == 'ordered' -> data is ordered_ranges list (each element is a ranges list)
      - mode == 'global' -> data is global_ranges list
    """
    lines = []
    with open(path, 'rt') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            lines.append(line)

    # If any line contains ':' treat file as per-id style
    if any(':' in ln for ln in lines):
        per_seq = {}
        for ln in lines:
            if ':' not in ln:
                # Non-per-id line in a per-id file: treat as global override? here we skip
                continue
            sid, rest = ln.split(':', 1)
            sid = sid.strip()
            rest = rest.strip()
            if rest == '':
                continue
            if rest == '-':
                per_seq[sid] = [None]
            else:
                per_seq[sid] = parse_positions_arg(rest)
        return 'per_id', per_seq

    # No ':' anywhere -> either ordered list or global single-line
    if len(lines) == 0:
        return 'global', []
    if len(lines) == 1:
        # single-line global ranges
        rngs = parse_positions_arg(lines[0])
        return 'global', rngs

    # multiple lines, no ':': treat as ordered per-sequence mapping
    ordered = []
    for ln in lines:
        if ln == '-':
            ordered.append([None])
        else:
            ordered.append(parse_positions_arg(ln))
    return 'ordered', ordered

def extract_segments(seq, ranges, tail_mode=False):
    """Return extracted subsequence for ranges. ranges may contain None meaning 'no cut'."""
    if not ranges:
        return ''
    if any(r is None for r in ranges):
        return seq
    L = len(seq)
    out_pieces = []
    for a, b in ranges:
        if tail_mode:
            start = max(1, b)
            end = L
        else:
            start = max(1, a)
            end = min(L, b)
        if start > end:
            continue
        out_pieces.append(seq[start-1:end])
    return ''.join(out_pieces)

def write_fasta_record(outfh, header, seq, wrap=0):
    outfh.write(f">{header}\n")
    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            outfh.write(seq[i:i+wrap] + "\n")
    else:
        outfh.write(seq + "\n")

def main():
    p = argparse.ArgumentParser(description="Extract subsequences from a FASTA using positions file.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--positions', help="Single positions string (global), e.g. '32-33' or '10-20,40-50'. Use '-' to mean no cut.")
    group.add_argument('--positions-file', help="File with positions. Either lines 'seq_id:32-33' or ordered list of ranges (one per FASTA record) or a single global line.")
    p.add_argument('--faa', required=True, help="Input FASTA (.faa)")
    p.add_argument('-o', '--out', default='positions_extracted.faa', help="Output FASTA")
    p.add_argument('--tail', action='store_true', help="Interpret each range A-B as 'take from B to end' (tail mode).")
    p.add_argument('--wrap', type=int, default=0, help="Wrap output sequences at given column (0 = no wrap).")
    p.add_argument('--verbose', action='store_true', help="Verbose messages")
    args = p.parse_args()

    if args.positions_file:
        mode, data = read_positions_file(args.positions_file)
    else:
        # single-line positions argument -> global
        mode = 'global'
        data = parse_positions_arg(args.positions)

    if args.verbose:
        print(f"Mode: {mode}", file=sys.stderr)
        if mode == 'global':
            print(f"Global ranges: {data}", file=sys.stderr)
        elif mode == 'per_id':
            print(f"Per-id ranges for {len(data)} ids", file=sys.stderr)
        else:
            print(f"Ordered ranges: {len(data)} lines", file=sys.stderr)

    written = 0
    seq_index = 0
    with open(args.out, 'wt') as outfh:
        for header, seq in fasta_iter(args.faa):
            seq_id = header.split()[0]
            seq_index += 1  # 1-based index for messages
            # determine ranges according to mode
            ranges = None
            if mode == 'per_id':
                ranges = data.get(seq_id, [])
            elif mode == 'global':
                ranges = data
            else:  # ordered
                idx = seq_index - 1
                if idx < len(data):
                    ranges = data[idx]
                else:
                    ranges = []  # no specification for this record

            if not ranges:
                if args.verbose:
                    print(f"Skipping {seq_id} (no ranges specified)", file=sys.stderr)
                continue

            subseq = extract_segments(seq, ranges, tail_mode=args.tail)
            if subseq == '':
                if args.verbose:
                    print(f"Warning: {seq_id} produced empty subsequence (ranges {ranges})", file=sys.stderr)
                continue
            write_fasta_record(outfh, header, subseq, wrap=args.wrap)
            written += 1

    if args.verbose:
        print(f"Done: wrote {written} sequences to {args.out}", file=sys.stderr)

if __name__ == '__main__':
    main()
"""
extract_by_positions.py

Trim/extract subsequences from a FASTA (.faa) file using position ranges.

Supports:
- Ordered per-line positions file (line N applies to FASTA record N). Use this when
  your pos.txt is a list of ranges (or '-' to mean "no cut") in the same order as
  the FASTA records.
- Per-sequence by-id entries like "seq_id:32-33" (overrides ordered/global).
- Global single-line positions applied to all sequences.

Coordinates are 1-based inclusive. --tail means A-B => take B..end.
"""
import argparse
import sys

def fasta_iter(path):
    """Yield (header_without_gt, sequence) for each FASTA record."""
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

def parse_range(r):
    parts = r.split('-', 1)
    if len(parts) != 2:
        raise ValueError(f"Bad range format: '{r}' (expected A-B)")
    a = int(parts[0])
    b = int(parts[1])
    return a, b

def parse_positions_arg(pos_arg):
    """Parse a positions string like '32-33' or '10-20,40-50' into list of (a,b) tuples or None markers."""
    pos_arg = pos_arg.strip()
    if not pos_arg:
        return []
    ranges = []
    for piece in pos_arg.split(','):
        piece = piece.strip()
        if not piece:
            continue
        if piece == '-':
            ranges.append(None)
            continue
        ranges.append(parse_range(piece))
    return ranges

def read_positions_file(path):
    """
    Read a positions file. Recognizes three forms:
      - Per-id: "seq_id:32-33" or "seq_id:10-20,40-50"  -> per_seq dict
      - Ordered list (no ':' anywhere): each line applies to successive FASTA records.
        A line '-' becomes [None] meaning "no cut".
      - Global single line: if file has a single (non-empty, non-comment) line without ':'
        it will be treated as the global ranges.
    Returns a tuple: (mode, data)
      - mode == 'per_id' -> data is per_seq dict (seq_id -> ranges list)
      - mode == 'ordered' -> data is ordered_ranges list (each element is a ranges list)
      - mode == 'global' -> data is global_ranges list
    """
    lines = []
    with open(path, 'rt') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            lines.append(line)

    # If any line contains ':' treat file as per-id style
    if any(':' in ln for ln in lines):
        per_seq = {}
        for ln in lines:
            if ':' not in ln:
                # Non-per-id line in a per-id file: treat as global override? here we skip
                continue
            sid, rest = ln.split(':', 1)
            sid = sid.strip()
            rest = rest.strip()
            if rest == '':
                continue
            if rest == '-':
                per_seq[sid] = [None]
            else:
                per_seq[sid] = parse_positions_arg(rest)
        return 'per_id', per_seq

    # No ':' anywhere -> either ordered list or global single-line
    if len(lines) == 0:
        return 'global', []
    if len(lines) == 1:
        # single-line global ranges
        rngs = parse_positions_arg(lines[0])
        return 'global', rngs

    # multiple lines, no ':': treat as ordered per-sequence mapping
    ordered = []
    for ln in lines:
        if ln == '-':
            ordered.append([None])
        else:
            ordered.append(parse_positions_arg(ln))
    return 'ordered', ordered

def extract_segments(seq, ranges, tail_mode=False):
    """Return extracted subsequence for ranges. ranges may contain None meaning 'no cut'."""
    if not ranges:
        return ''
    if any(r is None for r in ranges):
        return seq
    L = len(seq)
    out_pieces = []
    for a, b in ranges:
        if tail_mode:
            start = max(1, b)
            end = L
        else:
            start = max(1, a)
            end = min(L, b)
        if start > end:
            continue
        out_pieces.append(seq[start-1:end])
    return ''.join(out_pieces)

def write_fasta_record(outfh, header, seq, wrap=0):
    outfh.write(f">{header}\n")
    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            outfh.write(seq[i:i+wrap] + "\n")
    else:
        outfh.write(seq + "\n")

def main():
    p = argparse.ArgumentParser(description="Extract subsequences from a FASTA using positions file.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--positions', help="Single positions string (global), e.g. '32-33' or '10-20,40-50'. Use '-' to mean no cut.")
    group.add_argument('--positions-file', help="File with positions. Either lines 'seq_id:32-33' or ordered list of ranges (one per FASTA record) or a single global line.")
    p.add_argument('--faa', required=True, help="Input FASTA (.faa)")
    p.add_argument('-o', '--out', default='positions_extracted.faa', help="Output FASTA")
    p.add_argument('--tail', action='store_true', help="Interpret each range A-B as 'take from B to end' (tail mode).")
    p.add_argument('--wrap', type=int, default=0, help="Wrap output sequences at given column (0 = no wrap).")
    p.add_argument('--verbose', action='store_true', help="Verbose messages")
    args = p.parse_args()

    if args.positions_file:
        mode, data = read_positions_file(args.positions_file)
    else:
        # single-line positions argument -> global
        mode = 'global'
        data = parse_positions_arg(args.positions)

    if args.verbose:
        print(f"Mode: {mode}", file=sys.stderr)
        if mode == 'global':
            print(f"Global ranges: {data}", file=sys.stderr)
        elif mode == 'per_id':
            print(f"Per-id ranges for {len(data)} ids", file=sys.stderr)
        else:
            print(f"Ordered ranges: {len(data)} lines", file=sys.stderr)

    written = 0
    seq_index = 0
    with open(args.out, 'wt') as outfh:
        for header, seq in fasta_iter(args.faa):
            seq_id = header.split()[0]
            seq_index += 1  # 1-based index for messages
            # determine ranges according to mode
            ranges = None
            if mode == 'per_id':
                ranges = data.get(seq_id, [])
            elif mode == 'global':
                ranges = data
            else:  # ordered
                idx = seq_index - 1
                if idx < len(data):
                    ranges = data[idx]
                else:
                    ranges = []  # no specification for this record

            if not ranges:
                if args.verbose:
                    print(f"Skipping {seq_id} (no ranges specified)", file=sys.stderr)
                continue

            subseq = extract_segments(seq, ranges, tail_mode=args.tail)
            if subseq == '':
                if args.verbose:
                    print(f"Warning: {seq_id} produced empty subsequence (ranges {ranges})", file=sys.stderr)
                continue
            write_fasta_record(outfh, header, subseq, wrap=args.wrap)
            written += 1

    if args.verbose:
        print(f"Done: wrote {written} sequences to {args.out}", file=sys.stderr)

if __name__ == '__main__':
    main()
