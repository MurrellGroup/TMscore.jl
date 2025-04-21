module BioStructuresExt

using TMscore
using BioStructures

function TMscore.write_tempfile(arg::StructuralElementOrList, tempdir)
    path = joinpath(tempdir, "$(time_ns()).cif")
    writemmcif(path, arg)
    return path
end

end
