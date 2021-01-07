using Crystalline

nerrs = 0

for sgnum in 1:230
    lgirsd = get_lgirreps(sgnum, Val(3))
    has_printed_sgnum = false
    for (klab, lgirs) in lgirsd
        try
            characters.(lgirs)

        catch err
            if !has_printed_sgnum
                println(sgnum, " (", issymmorph(sgnum) ? "" : "non", "symmorphic; $(iuc(sgnum)))")
                has_printed_sgnum=true
            end
            print("   ", klab, ":   ")
            println(err.msg)
            global nerrs += 1
        end
    end
end

@info nerrs

## --------------------------------------------------------------------------------------- #
# Trying to find a point group

using Crystalline
using Crystalline: isapproxin

function try_find_negating_op(lg::LittleGroup, sg::SpaceGroup)
    kv = kvec(lg)
    ops = operations(lg)
    holosymmetric_pg = pointgroup("mmm") # TODO: hardcoded to sg 43 atm... 

    for op′ in unique([S"-x,-y,-z", sg..., holosymmetric_pg...])       
        kv′ = op′∘kv
        if kv′ ≈ -kv
            if all(op->isapproxin(inv(op′) ∘ op ∘ op′, lg), lg)
                return op′
            end
        else
            println("$(seitz(op′)) took $(string(kv)) to $(string(kv′))")
        end
    end
    return nothing
end

sgnum = 43 # first tricky case
lgs = collect(values(get_littlegroups(sgnum)))
sg  = spacegroup(sgnum)
try_find_negating_op.(lgs, Ref(sg))