using TMscore
using Test

using BioStructures
using ProteinChains

@testset "TMscore.jl" begin
    
    mktempdir() do dir
        downloadpdb(["1CRN", "1EJG"]; dir)
        downloadpdb(["1CRN", "1EJG"]; dir, format=MMCIFFormat)
        struc1 = read("$dir/1CRN.pdb", PDBFormat)
        struc2 = read("$dir/1EJG.pdb", PDBFormat)
        @test allequal([
            tmscore("$dir/1CRN.pdb", "$dir/1EJG.pdb"),
            tmscore("$dir/1CRN.cif", "$dir/1EJG.cif"),
            tmscore(struc1, struc2),
            tmscore(struc1, "$dir/1EJG.cif"),
            run_tmscore(struc1, struc2).tmscore,
            tmscore(pdb"1CRN", pdb"1EJG"),
        ])
    end

end
