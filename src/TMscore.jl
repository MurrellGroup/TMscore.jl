module TMscore

import TMscore_jll

export tmscore, run_tmscore, TMscoreResult

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

Base.show(io::IO, result::TMscoreResult) = print(io, TMscoreResult, "(\"\"\"\n", result.output, "\"\"\")")

"""
    TMscoreResult(output::String)::TMscoreResult

Parse the text printed by TMscore into a `TMscoreResult`.
"""
function TMscoreResult(output::String)::TMscoreResult
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
    run_tmscore(file1::String, file2::String; options...) -> TMscoreResult
    run_tmscore(struc1, struc2; options...) -> TMscoreResult

Invoke the TMscore binary and parse its output into a `TMscoreResult` struct.

`BioStructures.StructuralElementOrList` objects can also be used as input,
and will be written to a temporary file before processing.

## Options

### Boolean Flags
- `c::Bool`: Compare two complex structures with multiple chains (default: false)
- `seq::Bool`: Establish residue equivalence by sequence alignment instead of residue indices (default: false)
- `a::Bool`: TM-score normalized by the average length of two structures (default: false)
- `m::Bool`: Output TM-score rotation matrix (default: false)
- `fast::Bool`: Fast but slightly inaccurate alignment (default: false)
- `mirror::Bool`: Whether to align the mirror image of input structure (default: false)
- `het::Bool`: Whether to align residues marked as 'HETATM' in addition to 'ATOM' (default: false)

### Numeric Options
- `d::Number`: TM-score scaled by an assigned d0 (in Angstroms)
- `l::Int`: TM-score normalized by a specific length
- `ter::Int`: Strings to mark the end of a chain (default: 3)
  - 0: end of file
  - 1: ENDMDL or END
  - 2: ENDMDL, END, or different chain ID
  - 3: TER, ENDMDL, END or different chain ID
- `split::Int`: Whether to split PDB file into multiple chains (default: 0)
  - 0: treat the whole structure as one single chain
  - 1: treat each MODEL as a separate chain (-ter should be 0)
  - 2: treat each chain as a separate chain (-ter should be ≤1)
- `outfmt::Int`: Output format (default: 0)
  - 0: full output
  - 1: fasta format compact output
  - 2: tabular format very compact output
  - -1: full output, but without version or citation information
- `infmt1::Int`: Input format for chain1 (default: -1)
  - -1: automatically detect PDB or PDBx/mmCIF format
  - 0: PDB format
  - 1: SPICKER format
  - 2: xyz format
  - 3: PDBx/mmCIF format
- `infmt2::Int`: Input format for chain2 (default: -1, same options as infmt1)

### String Options
- `o::String`: Generate superposition output files with the given prefix
- `dir::String`: Perform all-against-all alignment among the list of PDB chains
- `dir1::String`: Use chain2 to search a list of PDB chains
- `dir2::String`: Use chain1 to search a list of PDB chains
- `suffix::String`: Add file name suffix to files listed by chain1_list or chain2_list (default: empty)
- `atom::String`: 4-character atom name used to represent a residue (default: " CA " for proteins, " C3'" for RNA/DNA)
- `mol::String`: Molecule type: RNA or protein (default: auto-detect)
"""
function run_tmscore(file1::AbstractString, file2::AbstractString; options...)
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

function run_tmscore(arg1, arg2; options...)
    mktempdir() do tempdir
        run_tmscore(write_tempfile(arg1, tempdir), write_tempfile(arg2, tempdir); options...)
    end
end

"""
    tmscore(file1::String, file2::String; options...) -> Float64
    tmscore(struc1, struc2; options...) -> Float64

Invoke the TMscore binary and return only the TM-score value.

This is a convenience function that calls `run_tmscore` internally and extracts
the `tmscore` field from the resulting `TMscoreResult`. See `run_tmscore` for
details on options and handling different input types.
"""
function tmscore(arg1, arg2; options...)
    result = run_tmscore(arg1, arg2; options...)
    return result.tmscore
end

end
