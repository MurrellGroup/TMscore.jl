# TMscore

[![Build Status](https://github.com/MurrellGroup/TMscore.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/TMscore.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installation

```julia
using Pkg
Pkg.add("TMscore")
```

## Usage

The `tmscore` function can take .pdb files, .cif files, BioStructures.jl objects, or any combination of these.

```julia
using TMscore
using BioStructures

downloadpdb(["1CRN", "1EJG"])
tmscore("1CRN.pdb", "1EJG.pdb") # 0.9975

downloadpdb(["1CRN", "1EJG"]; format=MMCIFFormat)
tmscore("1CRN.cif", "1EJG.cif") # 0.9975

struc1 = retrievepdb("1CRN"; format=MMCIFFormat)
struc2 = retrievepdb("1EJG"; format=MMCIFFormat)
tmscore(struc1, struc2) # 0.9975

tmscore(struc1, "1EJG.cif") # 0.9975
```

The `run_tmscore` function returns a `TMscoreResult` with pretty-printing for showing the original output from the binary:

```julia
julia> run_tmscore("1CRN.pdb", "1EJG.pdb")
TMscoreResult("""

 *************************************************************************
 *                                 TM-SCORE                              *
 * A scoring function to assess the similarity of protein structures     *
 * Based on statistics:                                                  *
 *       0.0 < TM-score < 0.17, random structural similarity             *
 *       0.5 < TM-score < 1.00, in about the same fold                   *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710 *
 * For comments, please email to: yangzhanglab@umich.edu                 *
 *************************************************************************

Structure1: 1CRN.pdb    Length=   46
Structure2: 1EJG.pdb    Length=   46 (by which all scores are normalized)
Number of residues in common=   46
RMSD of  the common residues=    0.106

TM-score    = 0.9975  (d0= 2.10)
MaxSub-score= 0.9991  (d0= 3.50)
GDT-TS-score= 1.0000 %(d<1)=1.0000 %(d<2)=1.0000 %(d<4)=1.0000 %(d<8)=1.0000
GDT-HA-score= 1.0000 %(d<0.5)=1.0000 %(d<1)=1.0000 %(d<2)=1.0000 %(d<4)=1.0000

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     -0.0392522747   0.9999885471  -0.0021262047   0.0042877688
 2     -0.0580775348   0.0021360752   0.9999950764  -0.0022987559
 3     -0.0385006946  -0.0042828601   0.0023078886   0.9999881653

Superposition in the TM-score: Length(d<5.0)= 46
(":" denotes the residue pairs of distance < 5.0 Angstrom)
TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN
::::::::::::::::::::::::::::::::::::::::::::::
TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN
1234567890123456789012345678901234567890123456

""")
```

The fields are as follows:

```julia
struct TMscoreResult
    output::String           # raw output from TMscore
    len1::Int                # length of structure 1
    len2::Int                # length of structure 2
    common::Int              # number of residues in common
    rmsd::Float64            # RMSD over the common residues
    tmscore::Float64         # TM‐score
    d0::Float64              # d0 distance
    maxsub::Float64          # MaxSub‐score
    maxsub_d0::Float64       # MaxSub d0
    gdt_ts::Float64          # GDT‐TS score
    gdt_ts_thresholds::Dict{Float64,Float64}  # (d<1, d<2, d<4, d<8) → value
    gdt_ha::Float64          # GDT‐HA score
    gdt_ha_thresholds::Dict{Float64,Float64}  # (d<0.5, d<1, d<2, d<4) → value
    rotation::Matrix{Float64}    # 3×3 rotation matrix
    translation::Vector{Float64} # length‐3 translation vector
end
```
