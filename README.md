# TMscore

[![Build Status](https://github.com/MurrellGroup/TMscore.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/TMscore.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Usage

The `tmscore` function can take .pdb files, .cif files, or even BioStructures.jl objects, or any combination of these:

```julia
using TMscore
using BioStructures

downloadpdb(["1CRN", "1EJG"])
tmscore("1CRN.pdb", "1EJG.pdb").tmscore # 0.9975

downloadpdb(["1CRN", "1EJG"]; format=MMCIFFormat)
tmscore("1CRN.cif", "1EJG.cif").tmscore # 0.9975

struc1 = retrievepdb("1CRN"; format=MMCIFFormat)
struc2 = retrievepdb("1EJG"; format=MMCIFFormat)
tmscore(struc1, struc2).tmscore # 0.9975

tmscore(struc1, "1EJG.cif").tmscore # 0.9975
```