using Crystalline, Test

@testset "Band representations" begin
allpaths = false
spinful  = false
# Check classification of spinless, time-reversal invariant band representations
# using Table 4 of Po, Watanabe, & Vishwanath Nature Commun 8, 50 (2017).
#    Z₂:          3, 11, 14, 27, 37, 48, 49, 50, 52, 53, 54, 56, 58, 60, 66, 
#                 68, 70, 75, 77, 82, 85, 86, 88, 103, 124, 128, 130, 162,
#                 163, 164, 165, 166, 167, 168, 171, 172, 176, 184, 192,
#                 201, 203
#    Z₂×Z₂:       12, 13, 15, 81, 84, 87
#    Z₂×Z₄:       147, 148
#    Z₂×Z₂×Z₂:    10, 83, 175
#    Z₂×Z₂×Z₂×Z₄: 2
@testset "Classification (TR invar., spinless)" begin
    table4 = ["Z₁" for _=1:230]
    table4[[3, 11, 14, 27, 37, 48, 49, 50, 52, 53, 54, 56, 58, 60, 66,
        68, 70, 75, 77, 82, 85, 86, 88, 103, 124, 128, 130, 162,
        163, 164, 165, 166, 167, 168, 171, 172, 176, 184, 192,
        201, 203]]                   .= "Z₂"
    table4[[12, 13, 15, 81, 84, 87]] .= "Z₂×Z₂"
    table4[[147, 148]]               .= "Z₂×Z₄"
    table4[[10, 83, 175]]            .= "Z₂×Z₂×Z₂"
    table4[2]                         = "Z₂×Z₂×Z₂×Z₄"

    for sgnum = 1:230
        BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
        @test classification(BRS) == table4[sgnum]
    end
end


# Check basis dimension of spinless, time-reversal invariant band representations
# using Table 2 of Po, Watanabe, & Vishwanath Nature Commun 8, 50 (2017).
@testset "Band structure dimen. (TR invar., spinless)" begin
    table4 = Vector{Int64}(undef,230)
    table4[[1, 4, 7, 9, 19, 29, 33, 76, 78, 144, 145, 169, 170]] .= 1
    table4[[8, 31, 36, 41, 43, 80, 92, 96, 110, 146, 161, 198]] .= 2
    table4[[5, 6, 18, 20, 26, 30, 32, 34, 40, 45, 46, 61, 106, 109, 151, 152, 153, 154, 159, 160, 171, 172, 173, 178, 179, 199, 212, 213]] .= 3
    table4[[24, 28, 37, 39, 60, 62, 77, 79, 91, 95, 102, 104, 143, 155, 157, 158, 185, 186, 196, 197, 210]] .= 4
    table4[[3, 14, 17, 27, 42, 44, 52, 56, 57, 94, 98, 100, 101, 108, 114, 122, 150, 156, 182, 214, 220]] .= 5
    table4[[11, 15, 35, 38, 54, 70, 73, 75, 88, 90, 103, 105, 107, 113, 142, 149, 167, 168, 184, 195, 205, 219]] .= 6
    table4[[13, 22, 23, 59, 64, 68, 82, 86, 117, 118, 120, 130, 163, 165, 180, 181, 203, 206, 208, 209, 211, 218, 228, 230]] .= 7
    table4[[21, 58, 63, 81, 85, 97, 116, 133, 135, 137, 148, 183, 190, 201, 217]] .= 8
    table4[[2, 25, 48, 50, 53, 55, 72, 99, 121, 126, 138, 141, 147, 188, 207, 216, 222]] .= 9
    table4[[12, 74, 93, 112, 119, 176, 177, 202, 204, 215]] .= 10
    table4[[66, 84, 128, 136, 166, 227]] .= 11
    table4[[51, 87, 89, 115, 129, 134, 162, 164, 174, 189, 193, 223, 226]] .= 12
    table4[[16, 67, 111, 125, 194, 224]] .= 13
    table4[[49, 140, 192, 200]] .= 14
    table4[[10, 69, 71, 124, 127, 132, 187]] .= 15
    table4[[225, 229]] .= 17
    table4[[65, 83, 131, 139, 175]] .= 18
    table4[221] = 22
    table4[191] = 24
    table4[[47, 123]] .= 27

    for sgnum = 1:230
        BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
        @test basisdim(BRS) == table4[sgnum]
    end
end
end