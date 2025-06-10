module BioStructuresExt

using TMscore
using BioStructures

function TMscore.write_tempfile(arg::MolecularStructure, tempdir)
    path = joinpath(tempdir, "$(arg.name)-$(rand(UInt)).cif")
    writemmcif(path, arg)
    return path
end

function TMscore.write_tempfile(arg::StructuralElementOrList, tempdir)
    path = joinpath(tempdir, "$(rand(UInt)).cif")
    writemmcif(path, arg)
    return path
end

end
