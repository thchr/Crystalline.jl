module RenderLatticeIsosurface
# -----------------------------------------------------------------------------------------
using Crystalline, PyPlot, MAT
using LinearAlgebra: norm
using Crystalline: °

export write_to_mat_file, write_to_mat_file, prepare_obj_files_with_matlab, 
    render_with_blender, render_lattice
# -----------------------------------------------------------------------------------------
# TODO: What about linux/mac?
const DEFAULT_BLENDERPATH = "C:\\Program Files\\Blender Foundation\\Blender\\blender.exe"
# -----------------------------------------------------------------------------------------

function write_to_mat_file(isoval, xyz, vals, Rs)
    matopen((@__DIR__)*"/isosurfvis/tmp/isosurfdata.mat", "w") do io
        write(io, "isoval", isoval)
        write(io, "xyz", collect(xyz))
        write(io, "vals", vals)
        write(io, "Rs", collect.(Rs))
    end
end
function prepare_obj_files_with_matlab()
    # create isosurfaces and isocaps from matlab; save to "isosurfaces.obj" and "isocaps.obj"
    folder = joinpath((@__DIR__), "isosurfvis")
    mfile = "calc_isosurface_export_obj"
    run(`cmd /c start /wait \"\" 
        matlab -wait -automation -noFigureWindows -sd $(folder) -r "$(mfile); exit"`, 
        wait=true)
end
function render_with_blender(render_savepath, render_engine="eevee")
    if lowercase(render_engine) ∉ ("eevee", "cycles", "blender_eevee")
        throw(DomainError(render_engine, "must be either 'eevee' or 'cycles'"))
    end

    blenderpath = get(ENV, "BLENDERPATH", DEFAULT_BLENDERPATH)
    pypath = joinpath((@__DIR__), "isosurfvis", "set_isosurface_obj.py")
    # The cmd command needs proper escaping of paths; a bit fragile - don't mess with it
    cmd = `cmd /c start /wait \"\"
           \"$(blenderpath)\" --background
           --python \"$(pypath)\"
           --
           \"$(render_savepath)\"
           $(render_engine)`
    run(Cmd(cmd, windows_verbatim=true), wait=true)
end

# clean up/remove /tmp/ directory
function clean_up()
    dir = joinpath((@__DIR__), "isosurfvis", "tmp")
    # remove contents of /tmp/ directory 
    foreach(p->rm(p), readdir(dir; join=true))
    # remove directory itself
    rm(dir; recursive=true, force=true)
end

"""
    render_lattice(flat, Rs::DirectBasis{3}, render_savepath::AbstractString; 
                   filling=nothing, isoval=nothing, render_engine::String="eevee", 
                   N::Integer=75, verbose::Bool=false)

Render a Fourier lattice `flat` with direct basis `Rs` using Blender, saving the output to
`render_savepath`.
"""
function render_lattice(flat, Rs::DirectBasis{3}, render_savepath::AbstractString; 
                        filling=nothing, isoval=nothing, render_engine::String="eevee", 
                        N::Integer=75, verbose::Bool=true, cleanup::Bool=true)

    filling === isoval === nothing && throw("Either filling or isoval must be given")
    (filling !== nothing && isoval !== nothing) && 
                            throw("filling and isoval cannot be supplied simultaneously")

    verbose && @info "Computing ($N×$N×$N) grid-discretized lattice data"
    if filling !== nothing
        xyz,vals,isoval = plot(flat, Rs; N=N, filling=filling)
    elseif isoval !== nothing
        xyz,vals,_ = plot(flat, Rs; N=N, isoval=isoval)
    end
    PyPlot.close("all") # TODO: Ought to not need to do this...

    if !isdir(joinpath(@__DIR__, "isosurfvis", "tmp"))  # ensure scratch directory exists
        mkpath(joinpath(@__DIR__, "isosurfvis", "tmp"))
    end 

    verbose && @info "Writing data to .mat file"
    write_to_mat_file(isoval, xyz, vals, Rs) # write .mat file w/ isosurface & Rs data

    verbose && @info "Generating .obj files for isosurface, isocaps, and unit cell in Matlab"
    prepare_obj_files_with_matlab()          # create .obj files to import in Blender

    verbose && @info "Rendering in Blender with $render_engine"
    render_with_blender(render_savepath, render_engine) # render in Blender

    verbose && @info "Done: cleaning up"
    cleanup && clean_up()
end

end # module