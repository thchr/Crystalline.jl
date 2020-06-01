using Crystalline, PyPlot, MAT, LinearAlgebra

# -----------------------------------------------------------------------------------------
# Functions

°(φ::Real) = deg2rad(φ)
function standard_directbasis(sgnum, D)
    D ≠ 3 && throw(DomainError(D, "Not implemented"))
    system = crystalsystem(sgnum, 3)

    if     system == "cubic"        # a=b=c & α=β=γ=90° (free: a)
        a = b = c = 1.0
        α = β = γ = °(90)
    elseif system == "hexagonal" || # a=b & α=β=90° & γ=120° (free: a,c)
           system == "trigonal"    
        a = b = 1.0;        c = .75
        α = β = °(90);      γ = °(120)
    elseif system == "tetragonal"   # a=b & α=β=γ=90° (free: a,c)
        a = b = 1.0;        c = .75
        α = β = γ = °(90)
    elseif system == "orthorhombic" # α=β=γ=90° (free: a,b,c)
        a = 1.0;            b, c = .75, .875
        α = β = γ = °(90)
    elseif system == "monoclinic"   # α=γ=90° (free: a,b,c,β≥90°)
        a = 1.0;            b, c = .75, .875
        α = γ = °(90);      β = °(100)
        if !Crystalline.isvalid_sphericaltriangle(α,β,γ)
            throw("Unexpected error")
        end
    elseif system == "triclinic"    # no conditions (free: a,b,c,α,β,γ)
        a = 1.0;            b, c = .75, .875
        α, β, γ = °(85), °(100), °(92.5)
        if !Crystalline.isvalid_sphericaltriangle(α,β,γ)
            throw("Unexpected error")
        end
    else 
        throw(DomainError(system))
    end        
    return Crystalline.crystal(a,b,c,α,β,γ)
end
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
function render_with_blender(render_savename, render_engine="eevee")
    if lowercase(render_engine) ∉ ("eevee", "cycles", "blender_eevee")
        throw(DomainError(render_engine, "must be either 'eevee' or 'cycles'"))
    end

    blenderpath = "C:\\Program Files\\Blender Foundation\\Blender\\blender.exe"
    pypath = joinpath((@__DIR__), "isosurfvis", "set_isosurface_obj.py")
    render_savepath = joinpath("renders", render_engine, render_savename)
    # The cmd command needs proper escaping of paths; a bit fragile - don't mess with it
    cmd = `cmd /c start /wait \"\"
           \"$(blenderpath)\" --background
           --python \"$(pypath)\"
           --
           $(render_savepath)
           $(render_engine)`
    run(Cmd(cmd, windows_verbatim=true), wait=true)
end

# clean up tmp/ directory
clean_up() = foreach(p->rm(p), readdir(joinpath((@__DIR__), "isosurfvis", "tmp")))


# -----------------------------------------------------------------------------------------
# Scripting

render  = true
#sgnum = 35; 
for sgnum in 1:230
    println(sgnum)
    dxmax = 2 .* (1,1,1)
    N     = 75              # number of points per dimension
    expon = 1.0             # exponent for coefficient fall-off
    filling = 0.5

    cntr  = centering(sgnum, 3)
    Rs    = standard_directbasis(sgnum, 3)
    flat  = levelsetlattice(sgnum, Val(3), idxmax)
    mflat = normscale!(modulate(flat), expon)
    #pflat = primitivize(mflat, cntr)
    #pRs   = primitivize(Rs, cntr)

    xyz,vals,isoval = plot(mflat, Rs; N=N, filling=filling)
    
    if render
        PyPlot.close("all")

        render_savename = "unitcell-sg$sgnum"
        render_engine   = "eevee"
        
        write_to_mat_file(isoval, xyz, vals, Rs) # write .mat file w/ isosurface & Rs data
        prepare_obj_files_with_matlab()          # create .obj files to import in Blender
        render_with_blender(render_savename, render_engine) # render in Blender
        #clean_up()
    end
end
