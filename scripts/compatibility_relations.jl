using Crystalline
using LightGraphs, MetaGraphs
# Add these packages to your active environment before running, if necessary
using Colors, GraphPlot, Compose

# --- utilities ---
function binary_layout(truefalses)
    posx = ifelse.(truefalses, -1.0, 1.0)
    posy = Vector{Float64}(undef, length(truefalses))
    Nᵗ = count(truefalses); Nᶠ = count(!, truefalses);
    if Nᵗ != 1
        posy[truefalses] = range(-1,1, length=Nᵗ)
    else 
        posy[truefalses] = [0.0]
    end
    if Nᶠ != 1
        posy[map(!, truefalses)] = range(-1,1, length=Nᶠ)
    else
        posy[map(!, truefalses)] = [0.0]
    end
    return posx, posy
end

# the `shell_layout` command in GraphPlot is very limited/bugged: fix it
function proper_shell_layout(G, nlist::Union{Nothing, Vector{Vector{Int}}} = nothing)
    GraphPlot._nv(G) == 1 && return [0.0], [0.0] # short circuit for one-element graph
    if nlist == nothing
        nlist = [collect(1:GraphPlot._nv(G)), Int[]]
    end
    radius = length(nlist[1]) > 1 ? 1.0 : 0.0
    locs_x = Vector{Float64}(undef, GraphPlot._nv(G))
    locs_y = Vector{Float64}(undef, GraphPlot._nv(G))
    for nodes in nlist
        θ = range(0, stop=2pi, length=length(nodes)+1)[1:end-1]
        locs_x[nodes] .= radius.*cos.(θ)
        locs_y[nodes] .= radius.*sin.(θ)
        radius += 1.0
    end   
    locs_x, locs_y
end

# --- graph visualization functions ---
function plot_special_to_nonspecial_connectivity(lgirvec) 
    sgnum = num(first(first(lgirvec)))
    D     = dim(first(first(lgirvec)))
    Nk    = length(lgirvec)

    _, kgraph = connectivity(lgirvec)

    special = Vector{Bool}(undef,Nk)
    cols = Vector{RGB}(undef, Nk)
    for idx in vertices(kgraph)
        kv = get_prop(kgraph, idx, :kvec)
        if isspecial(kv)  # special
            special[idx] = true
            cols[idx] = colorant"firebrick1"
        else              # non-special
            special[idx] = false
            cols[idx] = colorant"royalblue4"
        end
    end
    posx, posy = binary_layout(special)
    #posx, posy = proper_shell_layout(kgraph, [(1:Nk)[special], (1:Nk)[map(!, special)]])
    nodesize = [1+LightGraphs.outdegree(kgraph, v)+LightGraphs.indegree(kgraph, v) for v in LightGraphs.vertices(kgraph)]
    pg = gplot(kgraph, posx, posy,
               nodefillc=cols,
               nodelabel=" ".*get_prop.(Ref(kgraph), vertices(kgraph), Ref(:klab)).*" ".*
                         string.(get_prop.(Ref(kgraph), vertices(kgraph), Ref(:kvec))).*" ",
               edgestrokec=colorant"grey81",
               nodesize=nodesize,
               arrowlengthfrac=.025,
               nodelabeldist=1.0
               )
    draw(SVG("$(@__DIR__)/graphs/connectivity_special2nonspecial_$(D)d_sg$(sgnum).svg", 12cm, 12cm), pg)
    
    return nothing
end


# --- plotting the graph ---
function plot_connectivity(lgirvec) 
    sgnum = num(first(first(lgirvec)))
    D     = dim(first(first(lgirvec)))
    Nk    = length(lgirvec)

    cgraph, kgraph = Crystalline.connectivity(lgirvec)
    #posx, posy = proper_shell_layout(kgraph, [(1:Nk)[special], (1:Nk)[map(!, special)]])
    nodesize = [1+outdegree(cgraph, v)+indegree(cgraph, v) for v in collect(vertices(cgraph))]
    pg = gplot(cgraph, 
               nodefillc=colorant"royalblue4",
               nodelabel=get_prop.(Ref(cgraph), vertices(cgraph), Ref(:klab)).*"\n".*
                         string.(get_prop.(Ref(cgraph), vertices(cgraph), Ref(:kvec))), 
               edgelabel=join.(get_prop.(Ref(cgraph), edges(cgraph), Ref(:klabs)), ", "),
               edgestrokec=colorant"grey81",
               nodesize=nodesize,
               )
    draw(SVG("$(@__DIR__)/graphs/connectivity_$(D)d_sg$(sgnum).svg", 12cm, 12cm), pg)
    
    return nothing
end


# --- plot connectivity of special k-vectors ---
if false
    for sgnum in 1:MAX_SGNUM[3]
        plot_connectivity(get_lgirreps(sgnum, Val(3)))
    end
end