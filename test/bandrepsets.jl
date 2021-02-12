using Test
using Crystalline
import Crystalline: dlm2struct, normalizesubsup, formatirreplabel, subscriptify

Elementary_AllPaths_1="Wyckoff pos.|1a(1)|1a(1)
Band-Rep.|A↑G(1)|Aˢ↑G(1)
Decomposable|false|false
X:(1/2,0)|X₁(1)|Xˢ₂(1)
Y:(0,1/2)|Y₁(1)|Yˢ₂(1)
M:(1/2,1/2)|M₁(1)|Mˢ₂(1)
Γ:(0,0)|Γ₁(1)|Γˢ₂(1)
Ω:(u,v)|Ω₁(1)|Ωˢ₂(1)";  

Elementary_AllPaths_16="Wyckoff pos.|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|2b(3)|3c(2)|3c(2)|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|2b(3)|3c(2)|3c(2)
Band-Rep.|A↑G(1)|B↑G(1)|¹E₁↑G(1)|¹E₂↑G(1)|²E₁↑G(1)|²E₂↑G(1)|A₁↑G(2)|¹E↑G(2)|²E↑G(2)|A↑G(3)|B↑G(3)|¹Eˢ₁↑G(1)|¹Eˢ₂↑G(1)|¹Eˢ₃↑G(1)|²Eˢ₁↑G(1)|²Eˢ₂↑G(1)|²Eˢ₃↑G(1)|¹Eˢ↑G(2)|²Eˢ↑G(2)|Eˢ↑G(2)|¹Eˢ↑G(3)|²Eˢ↑G(3)
Decomposable|false|false|false|false|false|false|true|true|true|true|true|false|false|false|false|false|false|true|true|true|true|true
Γ:(0,0)|Γ₁(1)|Γ₂(1)|Γ₃(1)|Γ₄(1)|Γ₅(1)|Γ₆(1)|Γ₁(1)⊕Γ₂(1)|Γ₃(1)⊕Γ₄(1)|Γ₅(1)⊕Γ₆(1)|Γ₁(1)⊕Γ₃(1)⊕Γ₅(1)|Γ₂(1)⊕Γ₄(1)⊕Γ₆(1)|Γˢ₇(1)|Γˢ₉(1)|Γˢ₁₁(1)|Γˢ₈(1)|Γˢ₁₂(1)|Γˢ₁₀(1)|Γˢ₁₁(1)⊕Γˢ₁₂(1)|Γˢ₉(1)⊕Γˢ₁₀(1)|Γˢ₇(1)⊕Γˢ₈(1)|Γˢ₇(1)⊕Γˢ₉(1)⊕Γˢ₁₁(1)|Γˢ₈(1)⊕Γˢ₁₀(1)⊕Γˢ₁₂(1)
K:(1/3,1/3)|K₁(1)|K₁(1)|K₂(1)|K₂(1)|K₃(1)|K₃(1)|K₂(1)⊕K₃(1)|K₁(1)⊕K₃(1)|K₁(1)⊕K₂(1)|K₁(1)⊕K₂(1)⊕K₃(1)|K₁(1)⊕K₂(1)⊕K₃(1)|Kˢ₄(1)|Kˢ₅(1)|Kˢ₆(1)|Kˢ₄(1)|Kˢ₆(1)|Kˢ₅(1)|Kˢ₄(1)⊕Kˢ₅(1)|Kˢ₄(1)⊕Kˢ₆(1)|Kˢ₅(1)⊕Kˢ₆(1)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(1)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(1)
M:(1/2,0)|M₁(1)|M₂(1)|M₁(1)|M₂(1)|M₁(1)|M₂(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₂(1)|M₁(1)⊕2M₂(1)|2M₁(1)⊕M₂(1)|Mˢ₃(1)|Mˢ₃(1)|Mˢ₃(1)|Mˢ₄(1)|Mˢ₄(1)|Mˢ₄(1)|Mˢ₃(1)⊕Mˢ₄(1)|Mˢ₃(1)⊕Mˢ₄(1)|Mˢ₃(1)⊕Mˢ₄(1)|Mˢ₃(1)⊕2Mˢ₄(1)|2Mˢ₃(1)⊕Mˢ₄(1)
Ω:(u,v)|Ω₁(1)|Ω₁(1)|Ω₁(1)|Ω₁(1)|Ω₁(1)|Ω₁(1)|2Ω₁(1)|2Ω₁(1)|2Ω₁(1)|3Ω₁(1)|3Ω₁(1)|Ωˢ₂(1)|Ωˢ₂(1)|Ωˢ₂(1)|Ωˢ₂(1)|Ωˢ₂(1)|Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|3Ωˢ₂(1)|3Ωˢ₂(1)
Λ:(u,u)|Λ₁(1)|Λ₁(1)|Λ₁(1)|Λ₁(1)|Λ₁(1)|Λ₁(1)|2Λ₁(1)|2Λ₁(1)|2Λ₁(1)|3Λ₁(1)|3Λ₁(1)|Λˢ₂(1)|Λˢ₂(1)|Λˢ₂(1)|Λˢ₂(1)|Λˢ₂(1)|Λˢ₂(1)|2Λˢ₂(1)|2Λˢ₂(1)|2Λˢ₂(1)|3Λˢ₂(1)|3Λˢ₂(1)
Σ:(u,0)|Σ₁(1)|Σ₁(1)|Σ₁(1)|Σ₁(1)|Σ₁(1)|Σ₁(1)|2Σ₁(1)|2Σ₁(1)|2Σ₁(1)|3Σ₁(1)|3Σ₁(1)|Σˢ₂(1)|Σˢ₂(1)|Σˢ₂(1)|Σˢ₂(1)|Σˢ₂(1)|Σˢ₂(1)|2Σˢ₂(1)|2Σˢ₂(1)|2Σˢ₂(1)|3Σˢ₂(1)|3Σˢ₂(1)";


ElementaryTR_AllPaths_16="Wyckoff pos.|1a(6)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|3c(2)|3c(2)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|3c(2)
Band-Rep.|A↑G(1)|B↑G(1)|¹E₁²E₁↑G(2)|¹E₂²E₂↑G(2)|A₁↑G(2)|¹E²E↑G(4)|A↑G(3)|B↑G(3)|¹Eˢ₁²Eˢ₁↑G(2)|¹Eˢ₂²Eˢ₂↑G(2)|¹Eˢ₃²Eˢ₃↑G(2)|¹Eˢ²Eˢ↑G(4)|EˢEˢ↑G(4)|¹Eˢ²Eˢ↑G(6)
Decomposable|false|false|false|false|false|true|true|true|false|false|false|true|true|true
Γ:(0,0)|Γ₁(1)|Γ₂(1)|Γ₃Γ₅(2)|Γ₄Γ₆(2)|Γ₁(1)⊕Γ₂(1)|Γ₃Γ₅(2)⊕Γ₄Γ₆(2)|Γ₁(1)⊕Γ₃Γ₅(2)|Γ₂(1)⊕Γ₄Γ₆(2)|Γˢ₇Γˢ₈(2)|Γˢ₁₂Γˢ₉(2)|Γˢ₁₀Γˢ₁₁(2)|Γˢ₁₀Γˢ₁₁(2)⊕Γˢ₁₂Γˢ₉(2)|2Γˢ₇Γˢ₈(2)|Γˢ₇Γˢ₈(2)⊕Γˢ₁₀Γˢ₁₁(2)⊕Γˢ₁₂Γˢ₉(2)
K:(1/3,1/3)|K₁(1)|K₁(1)|K₂K₃(2)|K₂K₃(2)|K₂K₃(2)|2K₁(1)⊕K₂K₃(2)|K₁(1)⊕K₂K₃(2)|K₁(1)⊕K₂K₃(2)|2Kˢ₄(1)|Kˢ₅Kˢ₆(2)|Kˢ₅Kˢ₆(2)|2Kˢ₄(1)⊕Kˢ₅Kˢ₆(2)|2Kˢ₅Kˢ₆(2)|2Kˢ₄(1)⊕2Kˢ₅Kˢ₆(2)
M:(1/2,0)|M₁(1)|M₂(1)|2M₁(1)|2M₂(1)|M₁(1)⊕M₂(1)|2M₁(1)⊕2M₂(1)|M₁(1)⊕2M₂(1)|2M₁(1)⊕M₂(1)|Mˢ₃Mˢ₄(2)|Mˢ₃Mˢ₄(2)|Mˢ₃Mˢ₄(2)|2Mˢ₃Mˢ₄(2)|2Mˢ₃Mˢ₄(2)|3Mˢ₃Mˢ₄(2)
Ω:(u,v)|Ω₁(1)|Ω₁(1)|2Ω₁(1)|2Ω₁(1)|2Ω₁(1)|4Ω₁(1)|3Ω₁(1)|3Ω₁(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|4Ωˢ₂(1)|4Ωˢ₂(1)|6Ωˢ₂(1)
Λ:(u,u)|Λ₁(1)|Λ₁(1)|2Λ₁(1)|2Λ₁(1)|2Λ₁(1)|4Λ₁(1)|3Λ₁(1)|3Λ₁(1)|2Λˢ₂(1)|2Λˢ₂(1)|2Λˢ₂(1)|4Λˢ₂(1)|4Λˢ₂(1)|6Λˢ₂(1)
Σ:(u,0)|Σ₁(1)|Σ₁(1)|2Σ₁(1)|2Σ₁(1)|2Σ₁(1)|4Σ₁(1)|3Σ₁(1)|3Σ₁(1)|2Σˢ₂(1)|2Σˢ₂(1)|2Σˢ₂(1)|4Σˢ₂(1)|4Σˢ₂(1)|6Σˢ₂(1)";


Elementary_AllPaths_17="Wyckoff pos.|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|2b(3m)|3c(mm2)|3c(mm2)|3c(mm2)|3c(mm2)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|2b(3m)
Band-Rep.|A₁↑G(1)|A₂↑G(1)|B₁↑G(1)|B₂↑G(1)|E₁↑G(2)|E₂↑G(2)|A₁↑G(2)|A₂↑G(2)|E↑G(4)|A₁↑G(3)|A₂↑G(3)|B₁↑G(3)|B₂↑G(3)|Eˢ₁↑G(2)|Eˢ₂↑G(2)|Eˢ₃↑G(2)|¹Eˢ↑G(2)|²Eˢ↑G(2)|Eˢ₁↑G(4)
Decomposable|false|false|false|false|false|false|false|false|true|true|true|true|true|false|false|false|false|false|true
Γ:(0,0,0)|Γ₁(1)|Γ₂(1)|Γ₄(1)|Γ₃(1)|Γ₆(2)|Γ₅(2)|Γ₁(1)⊕Γ₄(1)|Γ₂(1)⊕Γ₃(1)|Γ₅(2)⊕Γ₆(2)|Γ₁(1)⊕Γ₅(2)|Γ₂(1)⊕Γ₅(2)|Γ₃(1)⊕Γ₆(2)|Γ₄(1)⊕Γ₆(2)|Γˢ₉(2)|Γˢ₈(2)|Γˢ₇(2)|Γˢ₇(2)|Γˢ₇(2)|Γˢ₈(2)⊕Γˢ₉(2)
K:(1/3,1/3)|K₁(1)|K₂(1)|K₂(1)|K₁(1)|K₃(2)|K₃(2)|K₃(2)|K₃(2)|K₁(1)⊕K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|Kˢ₆(2)|Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)|Kˢ₆(2)|Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(2)
M:(1/2,0)|M₁(1)|M₂(1)|M₄(1)|M₃(1)|M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₄(1)|M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₃(1)⊕M₄(1)|M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₄(1)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|2Mˢ₅(2)
Λ:(u,u)|Λ₁(1)|Λ₂(1)|Λ₂(1)|Λ₁(1)|Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕Λ₂(1)|2Λ₁(1)⊕2Λ₂(1)|2Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕2Λ₂(1)|2Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕2Λ₂(1)|Λˢ₃(1)⊕Λˢ₄(1)|Λˢ₃(1)⊕Λˢ₄(1)|Λˢ₃(1)⊕Λˢ₄(1)|Λˢ₃(1)⊕Λˢ₄(1)|Λˢ₃(1)⊕Λˢ₄(1)|2Λˢ₃(1)⊕2Λˢ₄(1)
Σ:(u,0)|Σ₁(1)|Σ₂(1)|Σ₁(1)|Σ₂(1)|Σ₁(1)⊕Σ₂(1)|Σ₁(1)⊕Σ₂(1)|2Σ₁(1)|2Σ₂(1)|2Σ₁(1)⊕2Σ₂(1)|2Σ₁(1)⊕Σ₂(1)|Σ₁(1)⊕2Σ₂(1)|Σ₁(1)⊕2Σ₂(1)|2Σ₁(1)⊕Σ₂(1)|Σˢ₃(1)⊕Σˢ₄(1)|Σˢ₃(1)⊕Σˢ₄(1)|Σˢ₃(1)⊕Σˢ₄(1)|Σˢ₃(1)⊕Σˢ₄(1)|Σˢ₃(1)⊕Σˢ₄(1)|2Σˢ₃(1)⊕2Σˢ₄(1)
Ω:(u,v)|Ω₁(1)|Ω₁(1)|Ω₁(1)|Ω₁(1)|2Ω₁(1)|2Ω₁(1)|2Ω₁(1)|2Ω₁(1)|4Ω₁(1)|3Ω₁(1)|3Ω₁(1)|3Ω₁(1)|3Ω₁(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|4Ωˢ₂(1)";

ElementaryTR_AllPaths_17="Wyckoff pos.|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|2b(3m)|3c(mm2)|3c(mm2)|3c(mm2)|3c(mm2)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|3c(mm2)
Band-Rep.|A₁↑G(1)|A₂↑G(1)|B₁↑G(1)|B₂↑G(1)|E₁↑G(2)|E₂↑G(2)|A₁↑G(2)|A₂↑G(2)|E↑G(4)|A₁↑G(3)|A₂↑G(3)|B₁↑G(3)|B₂↑G(3)|Eˢ₁↑G(2)|Eˢ₂↑G(2)|Eˢ₃↑G(2)|¹Eˢ²Eˢ↑G(4)|Eˢ₁↑G(4)|Eˢ↑G(6)
Decomposable|false|false|false|false|false|false|false|false|true|true|true|true|true|false|false|false|true|true|true
Γ:(0,0)|Γ₁(1)|Γ₂(1)|Γ₄(1)|Γ₃(1)|Γ₆(2)|Γ₅(2)|Γ₁(1)⊕Γ₄(1)|Γ₂(1)⊕Γ₃(1)|Γ₅(2)⊕Γ₆(2)|Γ₁(1)⊕Γ₅(2)|Γ₂(1)⊕Γ₅(2)|Γ₃(1)⊕Γ₆(2)|Γ₄(1)⊕Γ₆(2)|Γˢ₉(2)|Γˢ₈(2)|Γˢ₇(2)|2Γˢ₇(2)|Γˢ₈(2)⊕Γˢ₉(2)|Γˢ₇(2)⊕Γˢ₈(2)⊕Γˢ₉(2)
K:(1/3,1/3)|K₁(1)|K₂(1)|K₂(1)|K₁(1)|K₃(2)|K₃(2)|K₃(2)|K₃(2)|K₁(1)⊕K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|Kˢ₆(2)|Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)|2Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)⊕2Kˢ₆(2)
M:(1/2,0)|M₁(1)|M₂(1)|M₄(1)|M₃(1)|M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₄(1)|M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₃(1)⊕M₄(1)|M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₄(1)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|2Mˢ₅(2)|2Mˢ₅(2)|3Mˢ₅(2)
Λ:(u,u)|Λ₁(1)|Λ₂(1)|Λ₂(1)|Λ₁(1)|Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕Λ₂(1)|2Λ₁(1)⊕2Λ₂(1)|2Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕2Λ₂(1)|2Λ₁(1)⊕Λ₂(1)|Λ₁(1)⊕2Λ₂(1)|Λˢ₃(1)⊕Λˢ₄(1)|Λˢ₃(1)⊕Λˢ₄(1)|Λˢ₃(1)⊕Λˢ₄(1)|2Λˢ₃(1)⊕2Λˢ₄(1)|2Λˢ₃(1)⊕2Λˢ₄(1)|3Λˢ₃(1)⊕3Λˢ₄(1)
Σ:(u,0)|Σ₁(1)|Σ₂(1)|Σ₁(1)|Σ₂(1)|Σ₁(1)⊕Σ₂(1)|Σ₁(1)⊕Σ₂(1)|2Σ₁(1)|2Σ₂(1)|2Σ₁(1)⊕2Σ₂(1)|2Σ₁(1)⊕Σ₂(1)|Σ₁(1)⊕2Σ₂(1)|Σ₁(1)⊕2Σ₂(1)|2Σ₁(1)⊕Σ₂(1)|Σˢ₃(1)⊕Σˢ₄(1)|Σˢ₃(1)⊕Σˢ₄(1)|Σˢ₃(1)⊕Σˢ₄(1)|2Σˢ₃(1)⊕2Σˢ₄(1)|2Σˢ₃(1)⊕2Σˢ₄(1)|3Σˢ₃(1)⊕3Σˢ₄(1)
Ω:(u,v)|Ω₁(1)|Ω₁(1)|Ω₁(1)|Ω₁(1)|2Ω₁(1)|2Ω₁(1)|2Ω₁(1)|2Ω₁(1)|4Ω₁(1)|3Ω₁(1)|3Ω₁(1)|3Ω₁(1)|3Ω₁(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|2Ωˢ₂(1)|4Ωˢ₂(1)|4Ωˢ₂(1)|6Ωˢ₂(1)
";

Elementary_MaxPaths_16="Wyckoff pos.|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|2b(3)|3c(2)|3c(2)|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|2b(3)|3c(2)|3c(2)
Band-Rep.|A↑G(1)|B↑G(1)|¹E₁↑G(1)|¹E₂↑G(1)|²E₁↑G(1)|²E₂↑G(1)|A₁↑G(2)|¹E↑G(2)|²E↑G(2)|A↑G(3)|B↑G(3)|¹Eˢ₁↑G(1)|¹Eˢ₂↑G(1)|¹Eˢ₃↑G(1)|²Eˢ₁↑G(1)|²Eˢ₂↑G(1)|²Eˢ₃↑G(1)|¹Eˢ↑G(2)|²Eˢ↑G(2)|Eˢ↑G(2)|¹Eˢ↑G(3)|²Eˢ↑G(3)
Decomposable|false|false|false|false|false|false|true|true|true|true|true|false|false|false|false|false|false|true|true|true|true|true
Γ:(0,0)|Γ₁(1)|Γ₂(1)|Γ₃(1)|Γ₄(1)|Γ₅(1)|Γ₆(1)|Γ₁(1)⊕Γ₂(1)|Γ₃(1)⊕Γ₄(1)|Γ₅(1)⊕Γ₆(1)|Γ₁(1)⊕Γ₃(1)⊕Γ₅(1)|Γ₂(1)⊕Γ₄(1)⊕Γ₆(1)|Γˢ₇(1)|Γˢ₉(1)|Γˢ₁₁(1)|Γˢ₈(1)|Γˢ₁₂(1)|Γˢ₁₀(1)|Γˢ₁₁(1)⊕Γˢ₁₂(1)|Γˢ₉(1)⊕Γˢ₁₀(1)|Γˢ₇(1)⊕Γˢ₈(1)|Γˢ₇(1)⊕Γˢ₉(1)⊕Γˢ₁₁(1)|Γˢ₈(1)⊕Γˢ₁₀(1)⊕Γˢ₁₂(1)
K:(1/3,1/3)|K₁(1)|K₁(1)|K₂(1)|K₂(1)|K₃(1)|K₃(1)|K₂(1)⊕K₃(1)|K₁(1)⊕K₃(1)|K₁(1)⊕K₂(1)|K₁(1)⊕K₂(1)⊕K₃(1)|K₁(1)⊕K₂(1)⊕K₃(1)|Kˢ₄(1)|Kˢ₅(1)|Kˢ₆(1)|Kˢ₄(1)|Kˢ₆(1)|Kˢ₅(1)|Kˢ₄(1)⊕Kˢ₅(1)|Kˢ₄(1)⊕Kˢ₆(1)|Kˢ₅(1)⊕Kˢ₆(1)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(1)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(1)
M:(1/2,0)|M₁(1)|M₂(1)|M₁(1)|M₂(1)|M₁(1)|M₂(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₂(1)|M₁(1)⊕2M₂(1)|2M₁(1)⊕M₂(1)|Mˢ₃(1)|Mˢ₃(1)|Mˢ₃(1)|Mˢ₄(1)|Mˢ₄(1)|Mˢ₄(1)|Mˢ₃(1)⊕Mˢ₄(1)|Mˢ₃(1)⊕Mˢ₄(1)|Mˢ₃(1)⊕Mˢ₄(1)|Mˢ₃(1)⊕2Mˢ₄(1)|2Mˢ₃(1)⊕Mˢ₄(1)
";

Elementary_MaxPaths_17="Wyckoff pos.|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|2b(3m)|3c(mm2)|3c(mm2)|3c(mm2)|3c(mm2)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|2b(3m)
Band-Rep.|A₁↑G(1)|A₂↑G(1)|B₁↑G(1)|B₂↑G(1)|E₁↑G(2)|E₂↑G(2)|A₁↑G(2)|A₂↑G(2)|E↑G(4)|A₁↑G(3)|A₂↑G(3)|B₁↑G(3)|B₂↑G(3)|Eˢ₁↑G(2)|Eˢ₂↑G(2)|Eˢ₃↑G(2)|¹Eˢ↑G(2)|²Eˢ↑G(2)|Eˢ₁↑G(4)
Decomposable|false|false|false|false|false|false|false|false|true|true|true|true|true|false|false|false|false|false|true
Γ:(0,0,0)|Γ₁(1)|Γ₂(1)|Γ₄(1)|Γ₃(1)|Γ₆(2)|Γ₅(2)|Γ₁(1)⊕Γ₄(1)|Γ₂(1)⊕Γ₃(1)|Γ₅(2)⊕Γ₆(2)|Γ₁(1)⊕Γ₅(2)|Γ₂(1)⊕Γ₅(2)|Γ₃(1)⊕Γ₆(2)|Γ₄(1)⊕Γ₆(2)|Γˢ₉(2)|Γˢ₈(2)|Γˢ₇(2)|Γˢ₇(2)|Γˢ₇(2)|Γˢ₈(2)⊕Γˢ₉(2)
K:(1/3,1/3)|K₁(1)|K₂(1)|K₂(1)|K₁(1)|K₃(2)|K₃(2)|K₃(2)|K₃(2)|K₁(1)⊕K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|Kˢ₆(2)|Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)|Kˢ₆(2)|Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(2)
M:(1/2,0)|M₁(1)|M₂(1)|M₄(1)|M₃(1)|M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₄(1)|M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₃(1)⊕M₄(1)|M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₄(1)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|2Mˢ₅(2)
";

ElementaryTR_MaxPaths_16="Wyckoff pos.|1a(6)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|3c(2)|3c(2)|1a(6)|1a(6)|1a(6)|2b(3)|2b(3)|3c(2)
Band-Rep.|A↑G(1)|B↑G(1)|¹E₁²E₁↑G(2)|¹E₂²E₂↑G(2)|A₁↑G(2)|¹E²E↑G(4)|A↑G(3)|B↑G(3)|¹Eˢ₁²Eˢ₁↑G(2)|¹Eˢ₂²Eˢ₂↑G(2)|¹Eˢ₃²Eˢ₃↑G(2)|¹Eˢ²Eˢ↑G(4)|EˢEˢ↑G(4)|¹Eˢ²Eˢ↑G(6)
Decomposable|false|false|false|false|false|true|true|true|false|false|false|true|true|true
Γ:(0,0)|Γ₁(1)|Γ₂(1)|Γ₃Γ₅(2)|Γ₄Γ₆(2)|Γ₁(1)⊕Γ₂(1)|Γ₃Γ₅(2)⊕Γ₄Γ₆(2)|Γ₁(1)⊕Γ₃Γ₅(2)|Γ₂(1)⊕Γ₄Γ₆(2)|Γˢ₇Γˢ₈(2)|Γˢ₁₂Γˢ₉(2)|Γˢ₁₀Γˢ₁₁(2)|Γˢ₁₀Γˢ₁₁(2)⊕Γˢ₁₂Γˢ₉(2)|2Γˢ₇Γˢ₈(2)|Γˢ₇Γˢ₈(2)⊕Γˢ₁₀Γˢ₁₁(2)⊕Γˢ₁₂Γˢ₉(2)
K:(1/3,1/3)|K₁(1)|K₁(1)|K₂K₃(2)|K₂K₃(2)|K₂K₃(2)|2K₁(1)⊕K₂K₃(2)|K₁(1)⊕K₂K₃(2)|K₁(1)⊕K₂K₃(2)|2Kˢ₄(1)|Kˢ₅Kˢ₆(2)|Kˢ₅Kˢ₆(2)|2Kˢ₄(1)⊕Kˢ₅Kˢ₆(2)|2Kˢ₅Kˢ₆(2)|2Kˢ₄(1)⊕2Kˢ₅Kˢ₆(2)
M:(1/2,0)|M₁(1)|M₂(1)|2M₁(1)|2M₂(1)|M₁(1)⊕M₂(1)|2M₁(1)⊕2M₂(1)|M₁(1)⊕2M₂(1)|2M₁(1)⊕M₂(1)|Mˢ₃Mˢ₄(2)|Mˢ₃Mˢ₄(2)|Mˢ₃Mˢ₄(2)|2Mˢ₃Mˢ₄(2)|2Mˢ₃Mˢ₄(2)|3Mˢ₃Mˢ₄(2)
";

ElementaryTR_MaxPaths_17="Wyckoff pos.|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|2b(3m)|3c(mm2)|3c(mm2)|3c(mm2)|3c(mm2)|1a(6mm)|1a(6mm)|1a(6mm)|2b(3m)|2b(3m)|3c(mm2)
Band-Rep.|A₁↑G(1)|A₂↑G(1)|B₁↑G(1)|B₂↑G(1)|E₁↑G(2)|E₂↑G(2)|A₁↑G(2)|A₂↑G(2)|E↑G(4)|A₁↑G(3)|A₂↑G(3)|B₁↑G(3)|B₂↑G(3)|Eˢ₁↑G(2)|Eˢ₂↑G(2)|Eˢ₃↑G(2)|¹Eˢ²Eˢ↑G(4)|Eˢ₁↑G(4)|Eˢ↑G(6)
Decomposable|false|false|false|false|false|false|false|false|true|true|true|true|true|false|false|false|true|true|true
Γ:(0,0)|Γ₁(1)|Γ₂(1)|Γ₄(1)|Γ₃(1)|Γ₆(2)|Γ₅(2)|Γ₁(1)⊕Γ₄(1)|Γ₂(1)⊕Γ₃(1)|Γ₅(2)⊕Γ₆(2)|Γ₁(1)⊕Γ₅(2)|Γ₂(1)⊕Γ₅(2)|Γ₃(1)⊕Γ₆(2)|Γ₄(1)⊕Γ₆(2)|Γˢ₉(2)|Γˢ₈(2)|Γˢ₇(2)|2Γˢ₇(2)|Γˢ₈(2)⊕Γˢ₉(2)|Γˢ₇(2)⊕Γˢ₈(2)⊕Γˢ₉(2)
K:(1/3,1/3)|K₁(1)|K₂(1)|K₂(1)|K₁(1)|K₃(2)|K₃(2)|K₃(2)|K₃(2)|K₁(1)⊕K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|K₁(1)⊕K₃(2)|K₂(1)⊕K₃(2)|Kˢ₆(2)|Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)|2Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)⊕Kˢ₆(2)|Kˢ₄(1)⊕Kˢ₅(1)⊕2Kˢ₆(2)
M:(1/2,0)|M₁(1)|M₂(1)|M₄(1)|M₃(1)|M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)|M₁(1)⊕M₄(1)|M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₃(1)⊕M₄(1)|M₂(1)⊕M₃(1)⊕M₄(1)|M₁(1)⊕M₂(1)⊕M₃(1)|M₁(1)⊕M₂(1)⊕M₄(1)|Mˢ₅(2)|Mˢ₅(2)|Mˢ₅(2)|2Mˢ₅(2)|2Mˢ₅(2)|3Mˢ₅(2)
";

io_elementary_allpaths_16=IOBuffer(Elementary_AllPaths_16)
csv_elementary_allpaths_16=dlm2struct(io_elementary_allpaths_16, 16, true, false, false)

io_elementary_allpaths_17=IOBuffer(Elementary_AllPaths_17)
csv_elementary_allpaths_17=dlm2struct(io_elementary_allpaths_17, 17, true, false, false)

analytic_elementary_allpaths_16=calc_bandreps(16, Val(2), allpaths=true, timereversal=false)
analytic_elementary_allpaths_17=calc_bandreps(17, Val(2), allpaths=true, timereversal=false)

io_elementarytr_allpaths_16=IOBuffer(ElementaryTR_AllPaths_16)
csv_elementarytr_allpaths_16=dlm2struct(io_elementarytr_allpaths_16, 16, true, false, true)

io_elementarytr_allpaths_17=IOBuffer(ElementaryTR_AllPaths_17)
csv_elementarytr_allpaths_17=dlm2struct(io_elementarytr_allpaths_17, 17, true, false, true)

analytic_elementarytr_allpaths_16=calc_bandreps(16, Val(2), allpaths=true, timereversal=true)
analytic_elementarytr_allpaths_17=calc_bandreps(17, Val(2), allpaths=true, timereversal=true)

io_elementary_maxpaths_16=IOBuffer(Elementary_MaxPaths_16)
csv_elementary_maxpaths_16=dlm2struct(io_elementary_maxpaths_16, 16, false, false, false)

io_elementary_maxpaths_17=IOBuffer(Elementary_MaxPaths_17)
csv_elementary_maxpaths_17=dlm2struct(io_elementary_maxpaths_17, 17, false, false, false)

analytic_elementary_maxpaths_16=calc_bandreps(16, Val(2), allpaths=false, timereversal=false)
analytic_elementary_maxpaths_17=calc_bandreps(17, Val(2), allpaths=false, timereversal=false)

io_elementarytr_maxpaths_16=IOBuffer(ElementaryTR_MaxPaths_16)
csv_elementarytr_maxpaths_16=dlm2struct(io_elementarytr_maxpaths_16, 16, false, false, true)

io_elementarytr_maxpaths_17=IOBuffer(ElementaryTR_MaxPaths_17)
csv_elementarytr_maxpaths_17=dlm2struct(io_elementarytr_maxpaths_17, 17, false, false, true)

analytic_elementarytr_maxpaths_16=calc_bandreps(16, Val(2), allpaths=false, timereversal=true)
analytic_elementarytr_maxpaths_17=calc_bandreps(17, Val(2), allpaths=false, timereversal=true)


function wyckpos(br::BandRep)
    br.wyckpos
end

function irvec(br::BandRep)
    br.irvec
end

@testset "3D: Checking Wyckoff Position Sets" begin
    for sgnum in 1:87
        brs_a=calc_bandreps(sgnum, Val(3), allpaths=false, timereversal=true);
        brs_b=bandreps(sgnum, allpaths=false, timereversal=true);
        wyckpos_a=wyckpos.(brs_a)
        wyckpos_b=wyckpos.(brs_b)
        @test Set(wyckpos_a)==Set(wyckpos_b)
    end
end

@testset "3D: Checking dimension set" begin
    for sgnum in 1:87
        brs_a=calc_bandreps(sgnum, Val(3));
        brs_b=bandreps(sgnum);
        dims_a=dim.(brs_a)
        dims_b=dim.(brs_b)
        @test Set(dims_a)==Set(dims_b)
    end
end


@testset "3D: Checking dimension set without timereversal" begin
    for sgnum in 1:87
        brs_a=calc_bandreps(sgnum, Val(3), allpaths=false, timereversal=false);
        brs_b=bandreps(sgnum, allpaths=false, timereversal=false);
        dims_a=dim.(brs_a)
        dims_b=dim.(brs_b)
        @test Set(dims_a)==Set(dims_b)
    end
end

@testset "3D: Checking consistency of irlabs" begin
    for sgnum in 1:87
        brs_a=calc_bandreps(sgnum, Val(3), allpaths=false, timereversal=true);
        brs_b=bandreps(sgnum, allpaths=false, timereversal=true);
        irlabs_a=brs_a.irlabs
        irlabs_b=brs_b.irlabs
        @test Set(formatirreplabel.(irlabs_a))==Set(irlabs_b)
    end
end

@testset "3D: Checking dimensions of bandreps" begin
    for sgnum in 1:87
        brs_a=calc_bandreps(sgnum, Val(3), allpaths=false, timereversal=true);
        brs_b=bandreps(sgnum, allpaths=false, timereversal=true);
        dims_a=dim.(brs_a)
        dims_b=dim.(brs_b)
        @test Set(dims_a) == Set(dims_b)
    end
end

@testset "2D: Checking dims against known bandreps for allpaths=true, timereversal=false" begin

    dims_csv_elementary_allpaths_16=dim.(csv_elementary_allpaths_16)
    dims_analytic_elementary_allpaths_16=dim.(analytic_elementary_allpaths_16)

    dims_csv_elementary_allpaths_17=dim.(csv_elementary_allpaths_17)
    dims_analytic_elementary_allpaths_17=dim.(analytic_elementary_allpaths_17)

    @test Set(dims_csv_elementary_allpaths_16) ==  Set(dims_analytic_elementary_allpaths_16)
    @test Set(dims_csv_elementary_allpaths_17) ==  Set(dims_analytic_elementary_allpaths_17)

end

@testset "2D: Checking dims against known bandreps for allpaths=true, timereversal=true" begin
    
    dims_csv_elementarytr_allpaths_16=dim.(csv_elementarytr_allpaths_16)
    dims_analytic_elementarytr_allpaths_16=dim.(analytic_elementarytr_allpaths_16)

    dims_csv_elementarytr_allpaths_17=dim.(csv_elementarytr_allpaths_17)
    dims_analytic_elementarytr_allpaths_17=dim.(analytic_elementarytr_allpaths_17)

    @test Set(dims_csv_elementarytr_allpaths_16) ==  Set(dims_analytic_elementarytr_allpaths_16)
    @test Set(dims_csv_elementarytr_allpaths_17) ==  Set(dims_analytic_elementarytr_allpaths_17)

end

@testset "2D: Checking dims against known bandreps for allpaths=false, timereversal=false" begin

    dims_csv_elementary_maxpaths_16=dim.(csv_elementary_maxpaths_16)
    dims_analytic_elementary_maxpaths_16=dim.(analytic_elementary_maxpaths_16)

    dims_csv_elementary_maxpaths_17=dim.(csv_elementary_maxpaths_17)
    dims_analytic_elementary_maxpaths_17=dim.(analytic_elementary_maxpaths_17)

    @test Set(dims_csv_elementary_maxpaths_16) ==  Set(dims_analytic_elementary_maxpaths_16)
    @test Set(dims_csv_elementary_maxpaths_17) ==  Set(dims_analytic_elementary_maxpaths_17)

end

@testset "2D: dims Checking against known bandreps for allpaths=false, timereversal=true" begin
    
    dims_csv_elementarytr_maxpaths_16=dim.(csv_elementarytr_maxpaths_16)
    dims_analytic_elementarytr_maxpaths_16=dim.(analytic_elementarytr_maxpaths_16)

    dims_csv_elementarytr_maxpaths_17=dim.(csv_elementarytr_maxpaths_17)
    dims_analytic_elementarytr_maxpaths_17=dim.(analytic_elementarytr_maxpaths_17)

    @test Set(dims_csv_elementarytr_maxpaths_16) ==  Set(dims_analytic_elementarytr_maxpaths_16)
    @test Set(dims_csv_elementarytr_maxpaths_17) ==  Set(dims_analytic_elementarytr_maxpaths_17)

end


@testset "2D: Checking wyckoff positions against known bandreps for allpaths=true, timereversal=false" begin

    wyckpos_csv_elementary_allpaths_16=dim.(csv_elementary_allpaths_16)
    wyckpos_analytic_elementary_allpaths_16=dim.(analytic_elementary_allpaths_16)

    wyckpos_csv_elementary_allpaths_17=wyckpos.(csv_elementary_allpaths_17)
    wyckpos_analytic_elementary_allpaths_17=wyckpos.(analytic_elementary_allpaths_17)

    @test Set(wyckpos_csv_elementary_allpaths_16) ==  Set(wyckpos_analytic_elementary_allpaths_16)
    @test Set(wyckpos_csv_elementary_allpaths_17) ==  Set(wyckpos_analytic_elementary_allpaths_17)

end

@testset "2D: Checking consistency of irlabs" begin

    irlabs_csv_16=csv_elementary_maxpaths_16.irlabs
    irlabs_analytic_16=analytic_elementary_maxpaths_16.irlabs

    irlabs_csv_17=csv_elementary_maxpaths_17.irlabs
    irlabs_analytic_17=analytic_elementary_maxpaths_17.irlabs

    @test Set(formatirreplabel.(irlabs_csv_16))==Set(irlabs_analytic_16)
    @test Set(formatirreplabel.(irlabs_csv_17))==Set(irlabs_analytic_17)

end


@testset "2D: Checking irvecs against verified csv files" begin
    
    for irvec_analytic in irvec.(analytic_elementary_maxpaths_16)
        
        irvec_csv_idx=0

        for irvec_csv in irvec.(csv_elementary_maxpaths_16)
            if Set(irvec_csv)==Set(irvec_analytic)
                irvec_csv_idx+=1
            end
        end
        @test irvec_csv_idx !=0
    end

    for irvec_analytic in irvec.(analytic_elementary_maxpaths_17)
        
        irvec_csv_idx=0

        for irvec_csv in irvec.(csv_elementary_maxpaths_17)
            if Set(irvec_csv)==Set(irvec_analytic)
                irvec_csv_idx+=1
            end
        end
        @test irvec_csv_idx !=0
    end

end



@testset "2D: Checking irvecs against verified csv files- timereversal=true" begin
    
    for irvec_analytic in irvec.(analytic_elementarytr_maxpaths_16)
        
        irvec_csv_idx=0

        for irvec_csv in irvec.(csv_elementarytr_maxpaths_16)
            if Set(irvec_csv)==Set(irvec_analytic)
                irvec_csv_idx+=1
            end
        end
        @test irvec_csv_idx !=0
    end

    for irvec_analytic in irvec.(analytic_elementarytr_maxpaths_17)
        
        irvec_csv_idx=0

        for irvec_csv in irvec.(csv_elementarytr_maxpaths_17)
            if Set(irvec_csv)==Set(irvec_analytic)
                irvec_csv_idx+=1
            end
        end
        @test irvec_csv_idx !=0
    end

end


@testset "2D: Checking irvecs through irlab permutation" begin

    irlabs_csv_16=csv_elementary_maxpaths_16.irlabs
    irlabs_analytic_16=analytic_elementary_maxpaths_16.irlabs

    irvecs_analytic_16=irvec.(analytic_elementary_maxpaths_16)
    irvecs_csv_16=irvec.(csv_elementary_maxpaths_16)

    irvecs_analytic_17=irvec.(analytic_elementary_maxpaths_17)
    irvecs_csv_17=irvec.(csv_elementary_maxpaths_17)

    irlabs_csv_17=csv_elementary_maxpaths_17.irlabs
    irlabs_analytic_17=analytic_elementary_maxpaths_17.irlabs

    irlab_permuted_indices_16=Int64[] #Permuted indices of the irlabels (indices of the csv format for each analytic index in order)
    irlab_permuted_indices_17=Int64[] 

    for irlab in irlabs_analytic_16
        append!(irlab_permuted_indices_16, findfirst( ==(irlab), formatirreplabel.(irlabs_csv_16)))
    end

    for irlab in irlabs_analytic_17
        append!(irlab_permuted_indices_17, findfirst( ==(irlab), formatirreplabel.(irlabs_csv_17)))
    end
        
    for irvec_csv in irvecs_csv_16
        @test findfirst(==(irvec_csv[irlab_permuted_indices_16]), irvecs_analytic_16) != nothing
    end

    for irvec_csv in irvecs_csv_17
        @test findfirst(==(irvec_csv[irlab_permuted_indices_17]), irvecs_analytic_17) != nothing
    end

end


@testset "3D: Checking irvecs in 3D" begin

    for sgnum in 1:87

        brs_a=calc_bandreps(sgnum, Val(3), allpaths=false, timereversal=true);
        brs_b=bandreps(sgnum, allpaths=false, timereversal=true);

        irvecs_a = irvec.(brs_a)
        irvecs_b = irvec.(brs_b)

        for irvec_a in irvecs_a

            found_irvec=0

            for irvec_b in irvecs_b

                if Set(irvec_a)==Set(irvec_b)#findfirst(==(Set(irvec_a)), irvecs_b) != nothing
                    found_irvec +=1
                end
            end

            @test found_irvec != 0

        end

    end

end

@testset "Checking parent groups in 3D" begin

    parent_3d=[1, 3, 6, 7, 8, 25, 28, 32, 35, 75, 99, 100, 143, 156, 157, 168, 183]
    for sgnum in 1:17
        parent_group=parent_3d[sgnum]
        @test Set(dim.(calc_bandreps(sgnum, Val(2)))) ==  Set(dim.(bandreps(parent_group, 3)))
    end

end

@testset "Checking classification from parent groups in 3D" begin
    
    parent_3d=[1, 3, 6, 7, 8, 25, 28, 32, 35, 75, 99, 100, 143, 156, 157, 168, 183]
    for sgnum in 1:17
        parent_group=parent_3d[sgnum]
        @test classification(bandreps(parent_group)) == classification(calc_bandreps(sgnum, Val(2)))
    end

end




