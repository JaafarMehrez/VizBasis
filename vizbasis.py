#!/usr/bin/env python3
"""
Copyright 2026 Jaafar Mehrez
Email: jaafarmehrez@sjtu.edu.cn, jaafar@hpqc.org

Tolerant .gbs parser + radial overlay plots + numeric overlaps.
Usage:
    python vizbasis.py basisA.gbs basisB.gbs --element 
Outputs:
  ./basis_compare_plots/  -> PNG files (one per Angular Momentum, e.g. O_S.png)
  ./basis_compare_plots/compare_summary.txt -> numeric overlaps + diagnostics
"""
import re, os, argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

if hasattr(np, 'trapezoid'):
    trapz = np.trapezoid
else:
    trapz = np.trapz

RE_ELEMENT = re.compile(r'^\s*([A-Za-z]{1,2})\s+[-+]?\d+\s*$')
RE_SHELL  = re.compile(r'^\s*([A-Za-z]{1,2})\s+(\d+)\s+([+-]?\d+(\.\d+)?([eE][+-]?\d+)?)\s*$')

def read_text_normalized(path):
    b = Path(path).read_bytes()
    try:
        s = b.decode('utf-8')
    except:
        s = b.decode('latin-1')
    s = s.replace('\r', '\n')
    s = s.replace('\u00A0', ' ').replace('\u2009', ' ').replace('\u202F', ' ')
    s = s.replace('\t', ' ')
    return s

def parse_gbs_text(text, default_element=None):
    lines = text.splitlines()
    shells = []
    current_element = None
    i = 0
    N = len(lines)
    while i < N:
        raw = lines[i]
        i += 1
        if not raw.strip(): continue
        m_el = RE_ELEMENT.match(raw)
        if m_el:
            current_element = m_el.group(1)
            continue
        m_sh = RE_SHELL.match(raw)
        if m_sh:
            atype = m_sh.group(1).upper()
            nprim = int(m_sh.group(2))
            try: scale = float(m_sh.group(3))
            except: scale = 1.0
            primitives = []
            read = 0
            while read < nprim and i < N:
                rawp = lines[i].strip()
                i += 1
                if rawp == '': continue
                toks = rawp.split()
                try: nums = [float(t) for t in toks]
                except: i -= 1; break
                alpha = nums[0]
                coeffs = nums[1:] if len(nums) > 1 else [1.0]
                primitives.append((alpha, coeffs))
                read += 1
            element = current_element if current_element is not None else default_element
            shells.append({'element': element, 'atype': atype, 'nprim': nprim, 'scale': scale, 'primitives': primitives})
            continue
        parts = raw.split()
        if len(parts) >= 3 and parts[0].isalpha() and parts[1].isdigit():
            atype = parts[0].upper()
            nprim = int(parts[1])
            try: scale = float(parts[2])
            except: scale = 1.0
            primitives = []
            read = 0
            while read < nprim and i < N:
                rawp = lines[i].strip(); i += 1
                if rawp == '': continue
                toks = rawp.split()
                try: nums = [float(t) for t in toks]
                except: i -= 1; break
                alpha = nums[0]
                coeffs = nums[1:] if len(nums) > 1 else [1.0]
                primitives.append((alpha, coeffs)); read += 1
            element = current_element if current_element is not None else default_element
            shells.append({'element': element, 'atype': atype, 'nprim': nprim, 'scale': scale, 'primitives': primitives})
            continue
    return shells

L_MAP = {'S':0,'P':1,'D':2,'F':3,'G':4,'H':5,'I':6}
def build_contracted_radials(shell, r):
    # Map 'SP' to 0 (S) and 1 (P) logic is simplified here to just taking the first char
    # If using SP shells, usually GBS parsers treat them uniquely. 
    # For standard S, P, D blocks, first char is sufficient.
    l_char = shell['atype'][0]
    l = L_MAP.get(l_char, 0)
    
    if not shell['primitives']:
        return []
    
    ncols = max(len(p[1]) for p in shell['primitives'])
    out = []
    
    # Check if it's an SP shell (special case in some GBS formats)
    if shell['atype'] == 'SP' and ncols == 2:
        # Col 0 is S, Col 1 is P
        # S-part
        psi_s = np.zeros_like(r)
        for alpha, coeffs in shell['primitives']:
            eff_alpha = shell['scale'] * float(alpha)
            c = coeffs[0]
            psi_s += c * (r**0) * np.exp(-eff_alpha * r**2)
        out.append(('S', psi_s))
        # P-part
        psi_p = np.zeros_like(r)
        for alpha, coeffs in shell['primitives']:
            eff_alpha = shell['scale'] * float(alpha)
            c = coeffs[1]
            psi_p += c * (r**1) * np.exp(-eff_alpha * r**2)
        out.append(('P', psi_p))
    else:
        # Standard processing
        for col in range(ncols):
            psi = np.zeros_like(r)
            for alpha, coeffs in shell['primitives']:
                eff_alpha = shell['scale'] * float(alpha)
                c = coeffs[col] if col < len(coeffs) else 0.0
                psi += c * (r**l) * np.exp(-eff_alpha * r**2)
            out.append((l_char, psi))
    return out

def numeric_normalize(psi, r):
    integrand = 4*np.pi * (r**2) * (psi**2)
    norm = np.sqrt(trapz(integrand, r))
    if norm == 0:
        return psi, 0.0
    return psi / norm, norm

def radial_overlap(psi1, psi2, r):
    return trapz(4*np.pi * (r**2) * psi1 * psi2, r)

def vizbasis(fileA, fileB, element='O', outdir='basis_compare_plots', rmax=8.0, nr=3001):
    os.makedirs(outdir, exist_ok=True)
    textA = read_text_normalized(fileA)
    textB = read_text_normalized(fileB)
    shellsA = [s for s in parse_gbs_text(textA, default_element=element) if s['element'] == element]
    shellsB = [s for s in parse_gbs_text(textB, default_element=element) if s['element'] == element]

    summary = []
    summary.append(f"Parsed shells for element {element}: FileA={len(shellsA)} shells, FileB={len(shellsB)} shells\n")

    if len(shellsA)==0 or len(shellsB)==0:
        summary.append("No shells parsed in at least one file.")
        Path(outdir, "compare_summary.txt").write_text("\n".join(summary))
        print("Error: No shells found.")
        return

    r = np.linspace(1e-8, rmax, nr)
    
    dictA = defaultdict(list)
    dictB = defaultdict(list)

    def collect_radials(shell_list, target_dict):
        for s in shell_list:
            # build_contracted_radials, returns list of (L_label, psi_array)
            # to handle potential split shells (though rare in simple formats)
            results = build_contracted_radials(s, r)
            for (l_label, psi) in results:
                target_dict[l_label].append(psi)

    collect_radials(shellsA, dictA)
    collect_radials(shellsB, dictB)

    # Get all unique angular momenta present
    all_types = sorted(set(dictA.keys()) | set(dictB.keys()))

    for atype in all_types:
        radA = dictA.get(atype, [])
        radB = dictB.get(atype, [])
        
        # Normalize
        radA_n = [numeric_normalize(x, r)[0] for x in radA]
        radB_n = [numeric_normalize(x, r)[0] for x in radB]
        
        # Matrix calculations
        if radA_n and radB_n:
            ov = np.zeros((len(radA_n), len(radB_n)))
            l2 = np.zeros_like(ov)
            for i, a in enumerate(radA_n):
                for j, b in enumerate(radB_n):
                    ov[i,j] = radial_overlap(a, b, r)
                    diff = a - b
                    l2[i,j] = np.sqrt(trapz(4*np.pi * (r**2) * diff**2, r))
        else:
            ov = None
            l2 = None

        plt.style.use('seaborn-v0_8-paper')
        matplotlib.rcParams.update({
            'font.family': 'Georgia',
            'font.size': 12,
            'axes.labelsize': 13,
            'axes.titlesize': 15,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'legend.fontsize': 10,
            'figure.dpi': 300,
            'savefig.dpi': 300,
})
        plt.figure(figsize=(10, 5))

        # Subplot 1: Wavefunction
        plt.subplot(1, 2, 1)
        
        # Basis Set 1 (A) - Black Solid
        for i, a in enumerate(radA_n):
            lbl = "Basis Set 1" if i == 0 else None
            plt.plot(r, a, color='black', linestyle='-', label=lbl, linewidth=1.0)
            
        # Basis Set 2 (B) - Blue Dashed
        for j, b in enumerate(radB_n):
            lbl = "Basis Set 2" if j == 0 else None
            plt.plot(r, b, color='blue', linestyle='--', label=lbl, linewidth=1.0)
            
        plt.title(f"{element} {atype}-shell Orbitals")
        plt.xlabel("r (bohr)")
        plt.legend(loc='best', fontsize='small')
        plt.grid(False)

        # Subplot 2: Probability Density
        plt.subplot(1, 2, 2)
        for a in radA_n:
            plt.plot(r, 4*np.pi * r**2 * a**2, color='black', linestyle='-')
        for b in radB_n:
            plt.plot(r, 4*np.pi * r**2 * b**2, color='blue', linestyle='--')
            
        plt.title(f"{element} {atype} Radial Probability")
        plt.xlabel("r (bohr)")
        plt.grid(False)
        plt.tight_layout()

        fname = Path(outdir) / f"{element}_{atype}.png"
        plt.savefig(str(fname))
        plt.close()

        # Text Summary
        summary.append(f"\n=== Group {atype} ===")
        summary.append(f"Basis 1 count: {len(radA_n)}")
        summary.append(f"Basis 2 count: {len(radB_n)}")
        
        if ov is not None:
            # Check size to avoid massive text output
            if ov.size < 400: 
                summary.append("Overlap Matrix (Rows=Basis1, Cols=Basis2):")
                summary.append(np.array2string(np.round(ov, 4), separator=', '))
            else:
                summary.append("Overlap Matrix too large to display.")

    Path(outdir, "compare_summary.txt").write_text("\n".join(summary))
    print(f"Done. Processed groups: {all_types}. Results in '{outdir}'")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare two .gbs basis sets and visualize radial functions")
    p.add_argument('fileA'); p.add_argument('fileB')
    p.add_argument('--element','-e',default='O', help='element label (default O)')
    p.add_argument('--outdir',default='basis_compare_plots')
    p.add_argument('--rmax', type=float, default=8.0)
    args = p.parse_args()
    vizbasis(args.fileA, args.fileB, element=args.element, outdir=args.outdir, rmax=args.rmax)
