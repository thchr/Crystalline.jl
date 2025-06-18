using Crystalline
using Test

@testset "Symmetry eigenvalue analysis (`collect_compatible_symmetry_vectors`)" begin
    # slightly rounded symmetry eigenvalue data, corresponding to README ex. from MPBUtils
    symeigsv = [ # (plane group 10, photonic crystal data)
    [[1.0 + 0.0im, 0.99998 - 0.0im], [1.0 + 0.0im, -0.99998 + 0.0im], [1.0 + 0.0im, -0.99993 + 0.0im], [1.0 + 0.0im, -0.99995 + 0.0im], [1.0 + 0.0im, 0.99997 - 0.0im], [1.0 + 0.0im, 0.99978 - 0.0im], [1.0 + 0.0im, 0.99947 - 0.0im], [1.0 + 0.0im, 0.99962 - 0.0im], [1.0 + 0.0im, -0.99987 + 0.0im], [1.0 + 0.0im, 0.99995 - 0.0im]],
    [[1.0 + 0.0im, 0.99997 - 0.0im, 0.99998 - 0.0im, 0.99998 + 0.0im], [1.0 + 0.0im, -0.99992 + 0.0im, 1.0e-5 - 0.0im, 1.0e-5 + 0.0im], [1.0 + 0.0im, -0.99992 + 0.0im, -1.0e-5 + 0.0im, -1.0e-5 - 0.0im], [1.0 + 0.0im, 0.99998 - 0.0im, -0.99999 + 0.0im, -0.99999 - 0.0im], [1.0 + 0.0im, -0.99997 + 0.0im, -0.0 - 1.0e-5im, -0.0 + 1.0e-5im], [1.0 + 0.0im, -0.99997 + 0.0im, 0.0 + 1.0e-5im, 0.0 - 1.0e-5im], [1.0 + 0.0im, 0.99966 - 0.0im, -0.99982 - 0.0im, -0.99982 + 0.0im], [1.0 + 0.0im, 0.9996 - 0.0im, 0.99979 + 0.0im, 0.99979 - 0.0im], [1.0 + 0.0im, 0.99938 - 0.0im, -0.99969 + 0.0im, -0.99969 - 0.0im], [1.0 + 0.0im, -0.99921 + 0.0im, -0.00012 + 0.00593im, -0.00012 - 0.00593im]],
    [[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im], [1.0 + 0.0im, 0.99999 - 0.0im, 1.0 - 0.0im, 1.0 + 0.0im], [1.0 + 0.0im, -0.99996 + 0.0im, 1.0e-5 - 0.0im, 1.0e-5 + 0.0im], [1.0 + 0.0im, -0.99996 + 0.0im, -1.0e-5 + 0.0im, -1.0e-5 - 0.0im], [1.0 + 0.0im, 0.99985 + 0.0im, -0.99992 + 0.0im, -0.99992 - 0.0im], [1.0 + 0.0im, 0.99987 + 0.0im, -0.99993 + 0.0im, -0.99993 - 0.0im], [1.0 + 0.0im, 0.9997 + 0.0im, 0.99985 + 0.0im, 0.99985 - 0.0im], [1.0 + 0.0im, -0.99993 + 0.0im, 0.0 - 1.0e-5im, 0.0 + 1.0e-5im], [1.0 + 0.0im, -0.99993 - 0.0im, -0.0 + 1.0e-5im, -0.0 - 1.0e-5im], [1.0 + 0.0im, 0.99975 - 0.0im, -0.99988 - 0.0im, -0.99988 + 0.0im]]
    ]

    D, sgnum = 2, 10 # p4, with Z₂ indicator group
    brs = primitivize(calc_bandreps(sgnum, Val(D)))
    lgirsv = irreps(brs)

    ns_reference = parse.(SymmetryVector{2},
                          ["[X₁, M₁, Γ₁]",
                           "[3X₂, M₂+M₃M₄, Γ₁+Γ₃Γ₄]",
                           "[2X₁, M₃M₄, 2Γ₂]",
                           "[X₁, M₂, Γ₁]",
                           "[X₁+X₂, M₁+M₂, Γ₃Γ₄]"],
                          Ref(lgirsv))

    @test collect_compatible_symmetry_vectors(symeigsv, brs) == ns_reference
end # @testset
