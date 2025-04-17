module TMscore

import TMscore_jll

export tmscore, TMscoreResult

"""
    TMscoreResult

Holds all of the key metrics returned by the TMscore executable.
"""
struct TMscoreResult
    output::String             # raw output from TMscore
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

Base.string(result::TMscoreResult) = result.output
Base.show(io::IO, result::TMscoreResult) = print(io, TMscoreResult, "(\"\"\"\n", result.output, "\"\"\")")

"""
    parse_tmscore_output(output::String)::TMscoreResult

Parse the text printed by TMscore into a `TMscoreResult`.
"""
function TMscoreResult(output::String)
    len1 = 0; len2 = 0; common = 0; rmsd = 0.0
    tms = 0.0; d0 = 0.0
    ms = 0.0; ms_d0 = 0.0
    gdt_ts = 0.0; gdt_ha = 0.0
    gdt_ts_th = Dict{Float64,Float64}(); gdt_ha_th = Dict{Float64,Float64}()
    rotation = zeros(Float64, 3, 3)
    translation = zeros(Float64, 3)

    for line in split(output, '\n')
        s = strip(line)

        # Structure lengths
        if startswith(s, "Structure1:")
            m = match(r"Length=\s*(\d+)", s)
            if m !== nothing
                len1 = parse(Int, m.captures[1])
            end
        elseif startswith(s, "Structure2:")
            m = match(r"Length=\s*(\d+)", s)
            if m !== nothing
                len2 = parse(Int, m.captures[1])
            end

        # common residues & RMSD
        elseif startswith(s, "Number of residues in common=")
            m = match(r"common=\s*(\d+)", s)
            if m !== nothing
                common = parse(Int, m.captures[1])
            end
        elseif startswith(s, "RMSD")
            m = match(r"=\s*([\d\.]+)", s)
            if m !== nothing
                rmsd = parse(Float64, m.captures[1])
            end

        # TM-score + d0
        elseif startswith(s, "TM-score")
            m = match(r"TM-score\s*=\s*([\d\.]+).*\(d0=\s*([\d\.]+)\)", s)
            if m !== nothing
                tms = parse(Float64, m.captures[1])
                d0  = parse(Float64, m.captures[2])
            end

        # MaxSub-score + d0
        elseif startswith(s, "MaxSub-score")
            m = match(r"MaxSub-score\s*=\s*([\d\.]+).*\(d0=\s*([\d\.]+)\)", s)
            if m !== nothing
                ms   = parse(Float64, m.captures[1])
                ms_d0 = parse(Float64, m.captures[2])
            end

        # GDT-TS
        elseif startswith(s, "GDT-TS-score")
            m = match(r"=\s*([\d\.]+)", s)
            if m !== nothing
                gdt_ts = parse(Float64, m.captures[1])
            end
            for mm in eachmatch(r"\(d<([\d\.]+)\)=([\d\.]+)", s)
                thr = parse(Float64, mm.captures[1])
                val = parse(Float64, mm.captures[2])
                gdt_ts_th[thr] = val
            end

        # GDT-HA
        elseif startswith(s, "GDT-HA-score")
            m = match(r"=\s*([\d\.]+)", s)
            if m !== nothing
                gdt_ha = parse(Float64, m.captures[1])
            end
            for mm in eachmatch(r"\(d<([\d\.]+)\)=([\d\.]+)", s)
                thr = parse(Float64, mm.captures[1])
                val = parse(Float64, mm.captures[2])
                gdt_ha_th[thr] = val
            end

        end

        # rotation/translation rows
        m = match(r"^\s*(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", s)
        if m !== nothing
            i = parse(Int, m.captures[1])
            translation[i]     = parse(Float64, m.captures[2])
            rotation[i, 1:3]  .= parse.(Float64, m.captures[3:5])
        end
    end

    return TMscoreResult(
        output,
        len1, len2, common, rmsd,
        tms, d0, ms, ms_d0,
        gdt_ts, gdt_ts_th, gdt_ha, gdt_ha_th,
        rotation, translation
    )
end

tmscore_command(options::Vector{String}) = Cmd([TMscore_jll.TMscore().exec; options])

"""
    tmscore(file1::String, file2::String; options...)
    tmscore(struc1, struc2; options...)

Invoke the TMscore binary and parse its output.

`BioStructures.StructuralElementOrList` objects can also be used as input,
and will be written to a temporary file.

## Options

### Boolean Flags
- `c::Bool`: Compare two complex structures with multiple chains 
- `seq::Bool`: Establish residue equivalence by sequence alignment instead of residue indices
- `a::Bool`: TM-score normalized by the average length of two structures
- `m::Bool`: Output TM-score rotation matrix
- `fast::Bool`: Fast but slightly inaccurate alignment
- `mirror::Bool`: Whether to align the mirror image of input structure
- `het::Bool`: Whether to align residues marked as 'HETATM' in addition to 'ATOM'

### Numeric Options
- `d::Number`: TM-score scaled by an assigned d0 (in Angstroms)
- `l::Int`: TM-score normalized by a specific length
- `ter::Int`: Strings to mark the end of a chain
  - 0: No TER card
  - 1: TER separates different chains
  - 2: TER marks end of each chain
  - 3: EXIT separates different chains
- `split::Int`: Whether to split PDB file into multiple chains
  - 0: Don't split
  - 1: Split by chain ID
  - 2: Split by TER records
- `outfmt::Int`: Output format
  - 0: Full format
  - 1: Sequence and structure in fasta format
  - 2: Matrix format
  - -1: Compact format
- `infmt1::Int`: Input format for chain1
  - -1: Auto-detect
  - 0: PDB format
  - 1: SPICKER format
  - 2: xyz format
  - 3: FASTA format
- `infmt2::Int`: Input format for chain2 (same options as infmt1)

### String Options
- `o::String`: Generate superposition output files with the given prefix
- `dir::String`: Perform all-against-all alignment among the list of PDB chains
- `dir1::String`: Use chain2 to search a list of PDB chains
- `dir2::String`: Use chain1 to search a list of PDB chains
- `suffix::String`: Add file name suffix to files listed by chain1_list or chain2_list
- `atom::String`: 4-character atom name used to represent a residue
- `mol::String`: Molecule type: RNA or protein
"""
function tmscore(file1::AbstractString, file2::AbstractString; options...)
    cmd_vec = String[file1, file2]
    for (key, val) in pairs(options)
        if val === true
            push!(cmd_vec, "-$key")
        elseif val === false
            nothing
        else
            push!(cmd_vec, "-$key", string(val))
        end
    end

    cmd = tmscore_command(cmd_vec)
    stdout_buffer = IOBuffer()
    stderr_buffer = IOBuffer()
    success = try
        process = run(pipeline(cmd, stdout=stdout_buffer, stderr=stderr_buffer); wait=true)
        process.exitcode == 0
    catch e
        false
    end
    stdout_output = String(take!(stdout_buffer))
    stderr_output = String(take!(stderr_buffer))
    if success
        return TMscoreResult(stdout_output)
    else
        error_msg = "TMscore failed with error:\n" * 
                    (isempty(stderr_output) ? stdout_output : stderr_output)
        throw(ErrorException(error_msg))
    end
end

write_tempfile(filename::AbstractString, tempdir) = filename

function tmscore(x, y; options...)
    mktempdir() do tempdir
        tmscore(write_tempfile(x, tempdir), write_tempfile(y, tempdir); options...)
    end
end

end
