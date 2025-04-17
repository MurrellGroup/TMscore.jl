module BioStructuresExt

using TMscore
using BioStructures

function TMscore.write_tempfile(x::StructuralElementOrList, tempdir)
    path = joinpath(tempdir, "$(time_ns()).cif")
    writemmcif(path, x)
    return path
end

end
