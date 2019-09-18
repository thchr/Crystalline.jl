using SGOps, Test

#LGIRS=parselittlegroupirreps();

# see e.g. Bradley & Cracknell p. 151(bottom)-152(top)
@testset "Decomposition, order[star{k}]*order[G₀ᵏ] = order[G₀]" begin
for sgnum = 1:230
    cntr = centering(sgnum, 3)  
    # the "macroscopic" order is defined simply as the length of the 
    # point group associated with the space group
    sgops = operations(get_symops(sgnum, 3)) # from crawling Bilbao
    order_macroscopic = length(pointgroup(sgops))
    
    for kidx = 1:length(LGIRS[sgnum])
        kv = kvec(first(LGIRS[sgnum][kidx]))
        
        for iridx = 1:length(LGIRS[sgnum][kidx])
            lgir = LGIRS[sgnum][kidx][iridx]
            
            # number of k-vectors in the star of k
            order_starofk = length(starofk(sgops, kv, cntr)) 
            # number of operations in the little group of k
            order_pointgroupofk = length(pointgroup(operations(lgir)))

            # test that q⋅b=h, where
            #   q ≡ order[star{k}]
            #   b ≡ order[G₀ᵏ]      ≡ order of little co-group of k
            #   h ≡ order[G₀]       ≡ macroscopic order of G
            # and where G₀ denotes the point group derived from the space 
            # group G, Gᵏ denotes the little group of k derived from G, 
            # and G₀ᵏ denotes the point group derived from Gᵏ
            @test order_starofk*order_pointgroupofk == order_macroscopic
        end
    end
end
end

@testset "Macroscopic order, Bilbao vs. ISOTROPY" begin
for sgnum = 1:230
    sgops_bilbao = operations(get_symops(sgnum, 3))         # from crawling Bilbao
    sgops_isotropy = operations(first(first(LGIRS[sgnum]))) # from operations on Γ point irreps in ISOTROPY

    order_macroscopic_bilbao = length(pointgroup(sgops_bilbao))
    order_macroscopic_isotropy = length(pointgroup(sgops_isotropy))

    @test order_macroscopic_bilbao == order_macroscopic_isotropy
end   
end  
            