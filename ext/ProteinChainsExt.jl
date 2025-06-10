module ProteinChainsExt

using TMscore
using ProteinChains

function TMscore.write_tempfile(arg::Union{ProteinChain,ProteinStructure}, tempdir)
    path = joinpath(tempdir, "$(rand(UInt)).cif")
    writecif(path, arg)
    return path
end

end
