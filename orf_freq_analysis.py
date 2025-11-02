import math
from pathlib import Path
from collections import Counter
import argparse

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

DNA_STARTS = {"ATG", "GTG", "TTG"}
DNA_STOPS = {"TAA", "TAG", "TGA"}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
DIPEPTIDES = [a + b for a in AMINO_ACIDS for b in AMINO_ACIDS]

def revcomp(seq: str) -> str:
    tr = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(tr)[::-1]

def translate(nt: str) -> str:
    aa = []
    for i in range(0, len(nt) - 2, 3):
        cod = nt[i:i+3]
        a = codontab.get(cod, 'X')
        if a == '*': break
        aa.append(a)
    return ''.join(aa)

def read_fasta(path: Path):
    name, seqs = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name:
                    yield name, ''.join(seqs).upper().replace('U', 'T')
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if name:
            yield name, ''.join(seqs).upper().replace('U', 'T')

def find_orfs(seq: str, min_len: int):
    n = len(seq)
    orfs = []
    for frame in (0, 1, 2):
        starts, stops = [], []
        for i in range(frame, n - 2, 3):
            cod = seq[i:i+3]
            if cod in DNA_STARTS:
                starts.append(i)
            if cod in DNA_STOPS:
                stops.append(i)
        prev_stop = -3
        s_idx = 0
        for stop in stops:
            while s_idx < len(starts) and starts[s_idx] <= prev_stop:
                s_idx += 1
            if s_idx >= len(starts):
                break
            first = s_idx
            while s_idx < len(starts) and starts[s_idx] < stop:
                s_idx += 1
            if first < s_idx:
                start = starts[first]
                length = stop - start + 3
                if length >= min_len:
                    nt = seq[start:stop+3]
                    aa = translate(nt)
                    if aa:
                        orfs.append(aa)
            prev_stop = stop
    return orfs

def aa_freq(seqs):
    c = Counter()
    total = 0
    for s in seqs:
        for ch in s:
            if ch in AMINO_ACIDS:
                c[ch] += 1
                total += 1
    return {aa: c[aa] / total if total else 0 for aa in AMINO_ACIDS}

def dipep_freq(seqs):
    c = Counter()
    total = 0
    for s in seqs:
        for i in range(len(s) - 1):
            di = s[i:i+2]
            if di in DIPEPTIDES:
                c[di] += 1
                total += 1
    return {d: c[d] / total if total else 0 for d in DIPEPTIDES}

def cosine(a, b):
    dot = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na == 0 or nb == 0:
        return 1.0
    return 1 - dot / (na * nb)

def vec(keys, d):
    return [d.get(k, 0.0) for k in keys]

def write_phylip(names, D, path):
    with open(path, 'w') as f:
        f.write(f"{len(names)}\n")
        for i, n in enumerate(names):
            row = ' '.join(f"{x:.3f}" for x in D[i])
            f.write(f"{n}\t{row}\n")

def write_freqs(keys, rows, path):
    with open(path, 'w') as f:
        f.write('name\t' + '\t'.join(keys) + '\n')
        for name, data in rows.items():
            line = [f"{data.get(k, 0):.6f}" for k in keys]
            f.write(name + '\t' + '\t'.join(line) + '\n')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--glob', default='*.fasta')
    ap.add_argument('--min_orf_bp', type=int, default=100)
    ap.add_argument('--metric', default='cosine')
    ap.add_argument('--outdir', default='results')
    args = ap.parse_args()

    indir = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    all_freqs_aa = {}
    all_freqs_dp = {}
    names = []

    for f in sorted(indir.glob(args.glob)):
        for rec, seq in read_fasta(f):
            names.append(rec)
            all_orfs = find_orfs(seq, args.min_orf_bp) + find_orfs(revcomp(seq), args.min_orf_bp)
            all_freqs_aa[rec] = aa_freq(all_orfs)
            all_freqs_dp[rec] = dipep_freq(all_orfs)

    keys_aa = AMINO_ACIDS
    keys_dp = DIPEPTIDES

    D_aa = [[0]*len(names) for _ in names]
    D_dp = [[0]*len(names) for _ in names]
    for i in range(len(names)):
        for j in range(len(names)):
            if j < i:
                D_aa[i][j] = D_aa[j][i]
                D_dp[i][j] = D_dp[j][i]
            elif j == i:
                D_aa[i][j] = 0
                D_dp[i][j] = 0
            else:
                D_aa[i][j] = cosine(vec(keys_aa, all_freqs_aa[names[i]]),
                                    vec(keys_aa, all_freqs_aa[names[j]]))
                D_dp[i][j] = cosine(vec(keys_dp, all_freqs_dp[names[i]]),
                                    vec(keys_dp, all_freqs_dp[names[j]]))

    write_phylip(names, D_aa, outdir / 'distances_aa.phy')
    write_phylip(names, D_dp, outdir / 'distances_dipep.phy')
    write_freqs(keys_aa, all_freqs_aa, outdir / 'aa_freqs.tsv')
    write_freqs(keys_dp, all_freqs_dp, outdir / 'dipep_freqs.tsv')

    print('Analysis completed. Results can be found in the results folder: ', outdir)

if __name__ == '__main__':
    main()

