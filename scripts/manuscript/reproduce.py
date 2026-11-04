"""Regenerate five vector figures and validate the archived spectroscopy table."""
import argparse
from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/manuscript')
    args = parser.parse_args()
    output = args.output.resolve()
    commands = [
        ['check_seven_line_table.py'],
        ['make_figures.py', '--input', str(ROOT / 'benchmarks/data'),
         '--output', str(output)],
        ['make_level_schematic.py', '--output', str(output)],
    ]
    for script, *arguments in commands:
        subprocess.run([sys.executable, str(SCRIPT_DIR / script), *arguments], check=True)
    print(f'Five vector figures written to {output}')


if __name__ == '__main__':
    main()
