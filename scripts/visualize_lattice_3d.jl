using Pkg: activate, project
dirname(project().path) == (@__DIR__) || activate(@__DIR__)

using Crystalline

include((@__DIR__)*"/RenderLatticeIsosurface.jl")
using Main.RenderLatticeIsosurface
# -----------------------------------------------------------------------------------------
# Scripting

function default_directbasis(sgnum, D)
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

generate_savepath(savename, engine) = joinpath((@__DIR__), "isosurfvis", "renders", engine, savename)

render  = true
#sgnum = 35; 
for sgnum in 1:230
    println(sgnum)
    dxmax = 2 .* (1,1,1)
    N     = 75              # number of points per dimension
    expon = 1.0             # exponent for coefficient fall-off
    filling = 0.5

    cntr  = centering(sgnum, 3)
    Rs    = default_directbasis(sgnum, 3)
    flat  = levelsetlattice(sgnum, Val(3), idxmax)
    mflat = normscale!(modulate(flat), expon)
    #pflat = primitivize(mflat, cntr)
    #pRs   = primitivize(Rs, cntr)

    render_savename = "unitcell-sg$sgnum"
    render_engine   = "eevee"
    render_savepath = generate_savepath(render_savename, render_engine)
    render_lattice(mflat, Rs::DirectBasis{3}, render_savepath; 
                   filling=filling, render_engine=render_engine)
end
