using PrecompileTools: @setup_workload, @compile_workload
@setup_workload begin
    # __init__() doesn't run automatically during precompilation, so we need to do it 
    # invoke manually here; feels quite awful
    # (see also https://github.com/JuliaLang/PrecompileTools.jl/issues/32)
    __init__()

    @compile_workload begin
        function _precompile(sgnum, Dᵛ::Val{D}) where D
            sg = spacegroup(sgnum, Dᵛ)
            show(IOBuffer(), MIME"text/plain"(), sg)
            brs = calc_bandreps(sgnum, Dᵛ; timereversal=true, explicitly_real=true)
            show(IOBuffer(), MIME"text/plain"(), brs)
            lgirs = lgirreps(sgnum, Dᵛ)
            show(IOBuffer(), MIME"text/plain"(), lgirs)
        end
    _precompile(16, Val(2))
    _precompile(213, Val(3))
    end # @compile_workload
end # @setup_workload