"""Peptide feature encoders. Allele is fixed per task, so no MHC pseudoseq feature."""
import numpy as np

# BLOSUM62 matrix (20 standard AAs), standard values.
_AA = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62 = {
 'A':[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
 'R':[-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
 'N':[-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
 'D':[-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
 'C':[ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
 'Q':[-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
 'E':[-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
 'G':[ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
 'H':[-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
 'I':[-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
 'L':[-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
 'K':[-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
 'M':[-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
 'F':[-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
 'P':[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
 'S':[ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
 'T':[ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
 'W':[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
 'Y':[-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
 'V':[ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
}
_DEFAULT = [0]*20  # unknown AA (X) -> zeros

# Kyte-Doolittle hydropathy & charge for optional physchem features
_KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
 'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
 'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2,'X':0.0}
_CHARGE = {'D':-1,'E':-1,'K':1,'R':1,'H':0.1}

def blosum_encode(peptides, length=9):
    """BLOSUM62 row per residue, fixed length (9). -> (N, length*20) float array."""
    out = np.zeros((len(peptides), length*20), dtype=np.float32)
    for i, pep in enumerate(peptides):
        pep = pep[:length].ljust(length, 'X')  # 9-mers only in pass 1, but pad-safe
        for j, aa in enumerate(pep):
            out[i, j*20:(j+1)*20] = _BLOSUM62.get(aa, _DEFAULT)
    return out

def physchem(peptides):
    """length, mean hydropathy, net charge -> (N,3)."""
    out = np.zeros((len(peptides),3), dtype=np.float32)
    for i,pep in enumerate(peptides):
        out[i,0]=len(pep)
        out[i,1]=np.mean([_KD.get(a,0.0) for a in pep]) if pep else 0.0
        out[i,2]=sum(_CHARGE.get(a,0.0) for a in pep)
    return out
