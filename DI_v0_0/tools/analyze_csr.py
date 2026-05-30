#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict

def parse_file(path, max_top=20):
    dup_count = 0
    seen = set()
    zero_count = 0
    total = 0
    diag_vals = []
    top = []  # list of tuples (absval, row, col, val)

    with open(path, 'r', encoding='utf-8') as f:
        header = f.readline().strip()
        parts = header.split()
        try:
            nvars = int(parts[0])
            nnz_hdr = int(parts[1]) if len(parts) > 1 else None
        except Exception:
            # not a header, treat whole file as triples
            nvars = None
            nnz_hdr = None
            f.seek(0)

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                row = int(parts[0])
                col = int(parts[1])
                val = float(parts[2])
            except Exception:
                continue
            total += 1
            key = (row, col)
            if key in seen:
                dup_count += 1
            else:
                seen.add(key)
            if val == 0.0:
                zero_count += 1
            if row == col:
                diag_vals.append(val)
            a = abs(val)
            if len(top) < max_top:
                top.append((a, row, col, val))
                top.sort(reverse=True)
            else:
                if a > top[-1][0]:
                    top[-1] = (a, row, col, val)
                    top.sort(reverse=True)

    summary = {
        'nvars': nvars,
        'nnz_header': nnz_hdr,
        'entries_parsed': total,
        'unique_entries': len(seen),
        'duplicates_found': dup_count,
        'zero_coeffs': zero_count,
        'diag_count': len(diag_vals),
        'diag_min': min(diag_vals) if diag_vals else None,
        'diag_max': max(diag_vals) if diag_vals else None,
        'top': top,
    }
    return summary


def main():
    p = argparse.ArgumentParser(description='Analyze CSR triple file: row col val')
    p.add_argument('file', help='CSR file path (first line may be "nvars nnz")')
    args = p.parse_args()

    s = parse_file(args.file)
    print('CSR analysis summary:')
    print(f"  Header nvars={s['nvars']} nnz_header={s['nnz_header']}")
    print(f"  Parsed entries = {s['entries_parsed']}")
    print(f"  Unique (row,col) = {s['unique_entries']}")
    print(f"  Duplicate (repeated) entries = {s['duplicates_found']}")
    print(f"  Zero coefficients = {s['zero_coeffs']}")
    print(f"  Diagonal count = {s['diag_count']}")
    print(f"  Diagonal min/max = {s['diag_min']} / {s['diag_max']}")
    print('\nTop absolute entries:')
    for a,row,col,val in s['top']:
        print(f"  row={row} col={col} val={val} abs={a}")

if __name__ == '__main__':
    main()
