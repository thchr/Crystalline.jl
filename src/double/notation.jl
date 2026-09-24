# --- Mulliken notation for double-valued point group irreps ---

const PGIRLABS_CDML2MULLIKEN_3D_SPINFUL = Dict(
    # Mulliken labels of the double-valued irreps of the crystallographic point groups, keyed
    # by their CDML labels; the spinful counterpart of `PGIRLABS_CDML2MULLIKEN_3D`. From the
    # character tables of the double point groups of the Bilbao Crystallographic Server's
    # Representations DPG program (`representations_out.pl` with `tipogrupo=dbg`; see the
    # `dsg_crawl` artifact in `Artifacts.toml`), which list both labels. Bilbao's overline, marking a
    # double-valued irrep, is written as an appended `ˢ` (e.g., `²E₁gˢ`), as for the CDML
    # labels. Unlike `PGIRLABS_CDML2MULLIKEN_3D`, the labels are kept exactly as Bilbao gives
    # them (subscripts included, e.g. `E₁ˢ` in `312`), so that they agree with Bilbao's
    # spinful band representation tables.
    "1" => ImmutableDict("Γ₂ˢ"=>"Aˢ"),
    "-1" => ImmutableDict("Γ₂ˢ"=>"Aᵤˢ", "Γ₃ˢ"=>"Agˢ"),
    "2" => ImmutableDict("Γ₃ˢ"=>"²Eˢ", "Γ₄ˢ"=>"¹Eˢ"),
    "m" => ImmutableDict("Γ₃ˢ"=>"²Eˢ", "Γ₄ˢ"=>"¹Eˢ"),
    "2/m" => ImmutableDict("Γ₃ˢ"=>"²Egˢ", "Γ₄ˢ"=>"¹Egˢ", "Γ₅ˢ"=>"²Eᵤˢ", "Γ₆ˢ"=>"¹Eᵤˢ"),
    "222" => ImmutableDict("Γ₅ˢ"=>"Eˢ"),
    "mm2" => ImmutableDict("Γ₅ˢ"=>"Eˢ"),
    "mmm" => ImmutableDict("Γ₅ˢ"=>"Egˢ", "Γ₆ˢ"=>"Eᵤˢ"),
    "4" => ImmutableDict("Γ₅ˢ"=>"²E₂ˢ", "Γ₆ˢ"=>"²E₁ˢ", "Γ₇ˢ"=>"¹E₂ˢ", "Γ₈ˢ"=>"¹E₁ˢ"),
    "-4" => ImmutableDict("Γ₅ˢ"=>"²E₂ˢ", "Γ₆ˢ"=>"²E₁ˢ", "Γ₇ˢ"=>"¹E₂ˢ", "Γ₈ˢ"=>"¹E₁ˢ"),
    "4/m" => ImmutableDict("Γ₅ˢ"=>"²E₂gˢ", "Γ₆ˢ"=>"²E₁gˢ", "Γ₇ˢ"=>"¹E₂gˢ", "Γ₈ˢ"=>"¹E₁gˢ", "Γ₉ˢ"=>"²E₂ᵤˢ", "Γ₁₀ˢ"=>"²E₁ᵤˢ", "Γ₁₁ˢ"=>"¹E₂ᵤˢ", "Γ₁₂ˢ"=>"¹E₁ᵤˢ"),
    "422" => ImmutableDict("Γ₆ˢ"=>"E₂ˢ", "Γ₇ˢ"=>"E₁ˢ"),
    "4mm" => ImmutableDict("Γ₆ˢ"=>"E₂ˢ", "Γ₇ˢ"=>"E₁ˢ"),
    "-42m" => ImmutableDict("Γ₆ˢ"=>"E₂ˢ", "Γ₇ˢ"=>"E₁ˢ"),
    "-4m2" => ImmutableDict("Γ₆ˢ"=>"E₂ˢ", "Γ₇ˢ"=>"E₁ˢ"),
    "4/mmm" => ImmutableDict("Γ₆ˢ"=>"E₂gˢ", "Γ₇ˢ"=>"E₁gˢ", "Γ₈ˢ"=>"E₂ᵤˢ", "Γ₉ˢ"=>"E₁ᵤˢ"),
    "3" => ImmutableDict("Γ₄ˢ"=>"Eˢ", "Γ₅ˢ"=>"¹Eˢ", "Γ₆ˢ"=>"²Eˢ"),
    "-3" => ImmutableDict("Γ₄ˢ"=>"Egˢ", "Γ₅ˢ"=>"¹Egˢ", "Γ₆ˢ"=>"²Egˢ", "Γ₇ˢ"=>"Eᵤˢ", "Γ₈ˢ"=>"¹Eᵤˢ", "Γ₉ˢ"=>"²Eᵤˢ"),
    "312" => ImmutableDict("Γ₄ˢ"=>"²Eˢ", "Γ₅ˢ"=>"¹Eˢ", "Γ₆ˢ"=>"E₁ˢ"),
    "321" => ImmutableDict("Γ₄ˢ"=>"²Eˢ", "Γ₅ˢ"=>"¹Eˢ", "Γ₆ˢ"=>"E₁ˢ"),
    "3m1" => ImmutableDict("Γ₄ˢ"=>"²Eˢ", "Γ₅ˢ"=>"¹Eˢ", "Γ₆ˢ"=>"E₁ˢ"),
    "31m" => ImmutableDict("Γ₄ˢ"=>"²Eˢ", "Γ₅ˢ"=>"¹Eˢ", "Γ₆ˢ"=>"E₁ˢ"),
    "-31m" => ImmutableDict("Γ₄ˢ"=>"²Egˢ", "Γ₅ˢ"=>"¹Egˢ", "Γ₆ˢ"=>"²Eᵤˢ", "Γ₇ˢ"=>"¹Eᵤˢ", "Γ₈ˢ"=>"E₁gˢ", "Γ₉ˢ"=>"E₁ᵤˢ"),
    "-3m1" => ImmutableDict("Γ₄ˢ"=>"²Egˢ", "Γ₅ˢ"=>"¹Egˢ", "Γ₆ˢ"=>"²Eᵤˢ", "Γ₇ˢ"=>"¹Eᵤˢ", "Γ₈ˢ"=>"E₁gˢ", "Γ₉ˢ"=>"E₁ᵤˢ"),
    "6" => ImmutableDict("Γ₇ˢ"=>"²E₁ˢ", "Γ₈ˢ"=>"¹E₁ˢ", "Γ₉ˢ"=>"²E₂ˢ", "Γ₁₀ˢ"=>"¹E₃ˢ", "Γ₁₁ˢ"=>"²E₃ˢ", "Γ₁₂ˢ"=>"¹E₂ˢ"),
    "-6" => ImmutableDict("Γ₇ˢ"=>"²E₁ˢ", "Γ₈ˢ"=>"¹E₁ˢ", "Γ₉ˢ"=>"²E₂ˢ", "Γ₁₀ˢ"=>"¹E₃ˢ", "Γ₁₁ˢ"=>"²E₃ˢ", "Γ₁₂ˢ"=>"¹E₂ˢ"),
    "6/m" => ImmutableDict("Γ₇ˢ"=>"²E₁gˢ", "Γ₈ˢ"=>"¹E₁gˢ", "Γ₉ˢ"=>"²E₂gˢ", "Γ₁₀ˢ"=>"¹E₃gˢ", "Γ₁₁ˢ"=>"²E₃gˢ", "Γ₁₂ˢ"=>"¹E₂gˢ", "Γ₁₃ˢ"=>"²E₁ᵤˢ", "Γ₁₄ˢ"=>"¹E₁ᵤˢ", "Γ₁₅ˢ"=>"²E₂ᵤˢ", "Γ₁₆ˢ"=>"¹E₃ᵤˢ", "Γ₁₇ˢ"=>"²E₃ᵤˢ", "Γ₁₈ˢ"=>"¹E₂ᵤˢ"),
    "622" => ImmutableDict("Γ₇ˢ"=>"E₃ˢ", "Γ₈ˢ"=>"E₂ˢ", "Γ₉ˢ"=>"E₁ˢ"),
    "6mm" => ImmutableDict("Γ₇ˢ"=>"E₃ˢ", "Γ₈ˢ"=>"E₂ˢ", "Γ₉ˢ"=>"E₁ˢ"),
    "-62m" => ImmutableDict("Γ₇ˢ"=>"E₃ˢ", "Γ₈ˢ"=>"E₂ˢ", "Γ₉ˢ"=>"E₁ˢ"),
    "-6m2" => ImmutableDict("Γ₇ˢ"=>"E₃ˢ", "Γ₈ˢ"=>"E₂ˢ", "Γ₉ˢ"=>"E₁ˢ"),
    "6/mmm" => ImmutableDict("Γ₇ˢ"=>"E₃gˢ", "Γ₈ˢ"=>"E₂gˢ", "Γ₉ˢ"=>"E₁gˢ", "Γ₁₀ˢ"=>"E₃ᵤˢ", "Γ₁₁ˢ"=>"E₂ᵤˢ", "Γ₁₂ˢ"=>"E₁ᵤˢ"),
    "23" => ImmutableDict("Γ₅ˢ"=>"Eˢ", "Γ₆ˢ"=>"²Fˢ", "Γ₇ˢ"=>"¹Fˢ"),
    "m-3" => ImmutableDict("Γ₅ˢ"=>"Egˢ", "Γ₆ˢ"=>"²Fgˢ", "Γ₇ˢ"=>"¹Fgˢ", "Γ₈ˢ"=>"Eᵤˢ", "Γ₉ˢ"=>"²Fᵤˢ", "Γ₁₀ˢ"=>"¹Fᵤˢ"),
    "432" => ImmutableDict("Γ₆ˢ"=>"E₁ˢ", "Γ₇ˢ"=>"E₂ˢ", "Γ₈ˢ"=>"Fˢ"),
    "-43m" => ImmutableDict("Γ₆ˢ"=>"E₁ˢ", "Γ₇ˢ"=>"E₂ˢ", "Γ₈ˢ"=>"Fˢ"),
    "m-3m" => ImmutableDict("Γ₆ˢ"=>"E₁gˢ", "Γ₇ˢ"=>"E₂gˢ", "Γ₈ˢ"=>"E₁ᵤˢ", "Γ₉ˢ"=>"E₂ᵤˢ", "Γ₁₀ˢ"=>"Fgˢ", "Γ₁₁ˢ"=>"Fᵤˢ")
)

const PGIRLABS_CDML2MULLIKEN_3D_SPINFUL_COREP = Dict(
    # Same as `PGIRLABS_CDML2MULLIKEN_3D_SPINFUL`, but for the co-reps of the double-valued
    # irreps under time-reversal, as returned by `realify`; only point groups with such
    # co-reps are included, and only their co-reps (the remaining irreps are unchanged by
    # time-reversal). As for the single-valued co-reps (cf.
    # `PGIRLABS_CDML2MULLIKEN_3D_COREP`), a pair of complex conjugate irreps ¹Xˢ and ²Xˢ is
    # labelled Xˢ, unless Xˢ labels another irrep (as in `3` and `-3`); otherwise, a co-rep
    # is labelled by the concatenation of its irreps' labels (e.g., the doubled real irrep
    # Aˢ as AˢAˢ). Bilbao instead writes all co-reps as concatenations (e.g., `¹E₂ˢ²E₂ˢ`).
    "1" => ImmutableDict("Γ₂ˢΓ₂ˢ"=>"AˢAˢ"),
    "-1" => ImmutableDict("Γ₂ˢΓ₂ˢ"=>"AᵤˢAᵤˢ", "Γ₃ˢΓ₃ˢ"=>"AgˢAgˢ"),
    "2" => ImmutableDict("Γ₃ˢΓ₄ˢ"=>"Eˢ"),
    "m" => ImmutableDict("Γ₃ˢΓ₄ˢ"=>"Eˢ"),
    "2/m" => ImmutableDict("Γ₃ˢΓ₄ˢ"=>"Egˢ", "Γ₅ˢΓ₆ˢ"=>"Eᵤˢ"),
    "4" => ImmutableDict("Γ₅ˢΓ₇ˢ"=>"E₂ˢ", "Γ₆ˢΓ₈ˢ"=>"E₁ˢ"),
    "-4" => ImmutableDict("Γ₅ˢΓ₇ˢ"=>"E₂ˢ", "Γ₆ˢΓ₈ˢ"=>"E₁ˢ"),
    "4/m" => ImmutableDict("Γ₅ˢΓ₇ˢ"=>"E₂gˢ", "Γ₆ˢΓ₈ˢ"=>"E₁gˢ", "Γ₉ˢΓ₁₁ˢ"=>"E₂ᵤˢ", "Γ₁₀ˢΓ₁₂ˢ"=>"E₁ᵤˢ"),
    "3" => ImmutableDict("Γ₄ˢΓ₄ˢ"=>"EˢEˢ", "Γ₅ˢΓ₆ˢ"=>"¹Eˢ²Eˢ"),
    "-3" => ImmutableDict("Γ₄ˢΓ₄ˢ"=>"EgˢEgˢ", "Γ₅ˢΓ₆ˢ"=>"¹Egˢ²Egˢ", "Γ₇ˢΓ₇ˢ"=>"EᵤˢEᵤˢ", "Γ₈ˢΓ₉ˢ"=>"¹Eᵤˢ²Eᵤˢ"),
    "312" => ImmutableDict("Γ₄ˢΓ₅ˢ"=>"Eˢ"),
    "321" => ImmutableDict("Γ₄ˢΓ₅ˢ"=>"Eˢ"),
    "3m1" => ImmutableDict("Γ₄ˢΓ₅ˢ"=>"Eˢ"),
    "31m" => ImmutableDict("Γ₄ˢΓ₅ˢ"=>"Eˢ"),
    "-31m" => ImmutableDict("Γ₄ˢΓ₅ˢ"=>"Egˢ", "Γ₆ˢΓ₇ˢ"=>"Eᵤˢ"),
    "-3m1" => ImmutableDict("Γ₄ˢΓ₅ˢ"=>"Egˢ", "Γ₆ˢΓ₇ˢ"=>"Eᵤˢ"),
    "6" => ImmutableDict("Γ₇ˢΓ₈ˢ"=>"E₁ˢ", "Γ₉ˢΓ₁₂ˢ"=>"E₂ˢ", "Γ₁₀ˢΓ₁₁ˢ"=>"E₃ˢ"),
    "-6" => ImmutableDict("Γ₇ˢΓ₈ˢ"=>"E₁ˢ", "Γ₉ˢΓ₁₂ˢ"=>"E₂ˢ", "Γ₁₀ˢΓ₁₁ˢ"=>"E₃ˢ"),
    "6/m" => ImmutableDict("Γ₇ˢΓ₈ˢ"=>"E₁gˢ", "Γ₉ˢΓ₁₂ˢ"=>"E₂gˢ", "Γ₁₀ˢΓ₁₁ˢ"=>"E₃gˢ", "Γ₁₃ˢΓ₁₄ˢ"=>"E₁ᵤˢ", "Γ₁₅ˢΓ₁₈ˢ"=>"E₂ᵤˢ", "Γ₁₆ˢΓ₁₇ˢ"=>"E₃ᵤˢ"),
    "23" => ImmutableDict("Γ₆ˢΓ₇ˢ"=>"Fˢ"),
    "m-3" => ImmutableDict("Γ₆ˢΓ₇ˢ"=>"Fgˢ", "Γ₉ˢΓ₁₀ˢ"=>"Fᵤˢ")
)
