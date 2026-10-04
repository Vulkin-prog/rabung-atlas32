"""Explicitly compile the historical ternary catalogue worker; never run it."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nvcc', default='nvcc', help='Compatible NVIDIA CUDA C++ compiler')
    args = parser.parse_args()
    compiler = shutil.which(args.nvcc)
    if compiler is None:
        parser.exit(1, 'CUDA compiler unavailable; the stored-data checks need only Python.\n')
    source = root / 'src/ternary_catalogue'
    identity = json.loads((source / 'SOURCE_IDENTITY.json').read_text())
    output = root / 'build/ternary_catalogue/fast_catalogue'
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [compiler, *identity['historical_optimization'], '-std=c++17',
               *identity['historical_compile_options'], '-gencode', 'arch=compute_120,code=sm_120',
               str(source / 'fast_source/fast_worker.cu'), '-o', str(output)]
    subprocess.run(command, check=True)
    print('Compiled build/ternary_catalogue/fast_catalogue; no worker was started.')


if __name__ == '__main__':
    main()
