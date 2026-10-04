#!/usr/bin/env python3
"""Build source workers. CUDA compilation is opt-in; no worker is launched."""
import argparse
import json
from pathlib import Path
import subprocess
from atlas_io import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gpu', action='store_true', help='Additionally compile the historical CUDA SM120 worker')
    parser.add_argument('--build-dir', type=Path, default=ROOT / 'build')
    args = parser.parse_args(); args.build_dir.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'src/census'
    for name in ['small_cpu', 'prime_sieve', 'prime_count', 'moment_worker']:
        subprocess.run(['g++', '-O3', '-std=c++17', str(source / (name + '.cpp')),
                        '-o', str(args.build_dir / name)], check=True)
    if args.gpu:
        cfg = json.loads((source / 'GPU_CONFIG.json').read_text())['configuration']
        flags = ['nvcc', '-O3', '-std=c++17', '-lineinfo', '-gencode', 'arch=compute_120,code=sm_120']
        bridge = args.build_dir / 'baseline_bridge.o'
        subprocess.run(flags + ['-c', str(source / 'engine32/src/baseline_bridge.cu'), '-o', str(bridge)], check=True)
        subprocess.run(flags + cfg['compiler_options'] + ['-D' + v for v in cfg['macros']] +
                       [str(source / 'worker32.cu'), str(bridge), '-o', str(args.build_dir / 'worker32')], check=True)
    print('Workers built; no calculation launched.')


if __name__ == '__main__':
    main()
