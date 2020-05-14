# Conversion from the basis setting of ITA to ML/CDML: for an operator
# in ITA setting, gᴵᵀᴬ, the corresponding operator in ML/CDML setting,
# gᴹᴸ,  is obtained as gᴹᴸ = PgᴵᵀᴬP⁻¹.
#
# We can use this to transform an operator in CDML to one in the ITA 
# setting by gᴵᵀᴬ = P⁻¹gᴹᴸP; if P has rotation and rotation parts 𝐏 and 𝐩,
# we can achieve this by `gᴵᵀᴬ = transform(gᴹᴸ, P, p)`. [We use an opposite
# convention for P vs. P⁻¹ internally, so the switch from ML→ITA is correct.]
#
# The transformation matrices P are obtained from Stokes & Hatch's book, 
# Isotropy subgroups of the 230 crystallographic space groups, Table 6, 
# p. 6-1 to 6-4 (its use is described in their Section 4).
#
# If the settings agree between ITA and ML/CDML (i.e. if orientation and
# origo are the same; they are both in "conventional" settings, in the sense
# that they are not forced-primitive, as in B&C) the transformation is 
# trivial and we do not give one.

_unpack(pair::Pair{Int64,<:SymOperation}) = pair[1] => copy.(Crystalline.unpack(pair[2]))
const TRANSFORMS_CDML2ITA = Base.ImmutableDict(_unpack.([
    3:15   .=> Ref(SymOperation{3}("z,x,y"))
    38:41  .=> Ref(SymOperation{3}("-z,y,x"))
    48      => SymOperation{3}("x-1/4,y-1/4,z-1/4")
    50      => SymOperation{3}("x-1/4,y-1/4,z")
    59      => SymOperation{3}("x-1/4,y-1/4,z")
    68      => SymOperation{3}("x,y-1/4,z-1/4")
    70      => SymOperation{3}("x-7/8,y-7/8,z-7/8")
    85      => SymOperation{3}("x-3/4,y-1/4,z")
    86      => SymOperation{3}("x-3/4,y-3/4,z-3/4")
    88      => SymOperation{3}("x,y-3/4,z-7/8")
    125     => SymOperation{3}("x-3/4,y-3/4,z")
    126     => SymOperation{3}("x-3/4,y-3/4,z-3/4")
    129     => SymOperation{3}("x-3/4,y-1/4,z")
    130     => SymOperation{3}("x-3/4,y-1/4,z")
    133     => SymOperation{3}("x-3/4,y-1/4,z-3/4")
    134     => SymOperation{3}("x-3/4,y-1/4,z-3/4")
    137     => SymOperation{3}("x-3/4,y-1/4,z-3/4")
    138     => SymOperation{3}("x-3/4,y-1/4,z-3/4")
    141     => SymOperation{3}("x,y-1/4,z-7/8")
    142     => SymOperation{3}("x,y-1/4,z-7/8")
    146     => SymOperation{3}("y,-x+y,z")
    148     => SymOperation{3}("y,-x+y,z")
    151     => SymOperation{3}("x,y,z-1/6")
    153     => SymOperation{3}("x,y,z-5/6")
    154     => SymOperation{3}("x,y,z-1/6")
    155     => SymOperation{3}("y,-x+y,z")
    160     => SymOperation{3}("y,-x+y,z")
    161     => SymOperation{3}("y,-x+y,z")
    166     => SymOperation{3}("y,-x+y,z")
    167     => SymOperation{3}("y,-x+y,z")
    201     => SymOperation{3}("x-3/4,y-3/4,z-3/4")
    203     => SymOperation{3}("x-7/8,y-7/8,z-7/8")
    222     => SymOperation{3}("x-3/4,y-3/4,z-3/4")
    224     => SymOperation{3}("x-3/4,y-3/4,z-3/4")
])...)