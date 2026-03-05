# VizBasis

A Python tool for comparing Gaussian basis set (.gbs) files through radial overlay plots and numeric overlap calculations.

## Features

- Tolerant `.gbs` parser
- Radial wavefunction visualization
- Radial probability density plots
- Numeric overlap matrix computation
- Handles SP shells (combined S and P)

## Installation

Requires Python 3 with the following dependencies:

```bash
pip install numpy matplotlib
```

## Usage

```bash
python vizbasis.py basisA.gbs basisB.gbs --element O
```

### Arguments

- `fileA` - First basis set file (.gbs)
- `fileB` - Second basis set file (.gbs)
- `--element`, `-e` - Element symbol to analyze (default: O)
- `--outdir` - Output directory (default: basis_compare_plots)
- `--rmax` - Maximum radial distance in bohr (default: 8.0)

### Output

- `./basis_compare_plots/` - PNG files (one per angular momentum, e.g., `O_S.png`, `O_P.png`)
- `./basis_compare_plots/compare_summary.txt` - Numeric overlaps and diagnostics

## Example

```bash
python vizbasis.py example/basis-A example/basis-B --element O
```

## License

GNU General Public License v3.0 - See LICENSE file for details.

## Author

Jaafar Mehrez  
Email: jaafarmehrez@sjtu.edu.cn, jaafar@hpqc.org
