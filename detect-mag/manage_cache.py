#!/usr/bin/env python3
"""List cache usage; optionally remove known oversized intermediates.

Stop all notebooks using this directory before using --delete. This tool cannot
determine whether a remote kernel still holds a file open.
"""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', type=Path, default=Path(__file__).resolve().parent / 'cache')
    parser.add_argument('--threshold-gib', type=float, default=1.0)
    parser.add_argument('--delete', action='store_true', help='Delete listed known intermediates; kernels must be stopped')
    args = parser.parse_args()
    if args.threshold_gib <= 0:
        parser.error('--threshold-gib must be positive')
    if not args.cache_dir.is_dir():
        parser.error(f'Cache directory does not exist: {args.cache_dir}')
    files = sorted((p for p in args.cache_dir.iterdir() if p.is_file() and not p.is_symlink()),
                   key=lambda p: p.stat().st_size, reverse=True)
    total = sum(p.stat().st_size for p in files)
    reclaimed = 0
    for path in files:
        size = path.stat().st_size
        known = ((path.name.startswith('theory_rm_grid_') and path.suffix == '.npy')
                 or (path.name.startswith(('coeval_work_', 'theory_work_'))
                     and path.name.endswith('.tmp.npy')))
        candidate = known and size > args.threshold_gib * 2**30
        action = 'DELETE' if args.delete and candidate else 'candidate' if candidate else 'keep'
        print(f'{size / 2**30:8.3f} GiB  {action:9s} {path.name}')
        if args.delete and candidate:
            path.unlink()
            reclaimed += size
    print(f'Total before: {total / 2**30:.3f} GiB; removed: {reclaimed / 2**30:.3f} GiB')
    if not args.delete:
        print('Read-only listing. Stop notebook kernels, then use --delete to remove candidates.')


if __name__ == '__main__':
    main()
