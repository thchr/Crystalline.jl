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

_unpack(pair::Pair{Int64,SymOperation}) = pair[1] => copy.(SGOps.unpack(pair[2]))
const TRANSFORMS_CDML2ITA = Base.ImmutableDict(_unpack.([
    3:15   .=> Ref(SymOperation("z,x,y"))
    38:41  .=> Ref(SymOperation("-z,y,x"))
    48      => SymOperation("x-1/4,y-1/4,z-1/4")
    50      => SymOperation("x-1/4,y-1/4,z")
    59      => SymOperation("x-1/4,y-1/4,z")
    68      => SymOperation("x,y-1/4,z-1/4")
    70      => SymOperation("x-7/8,y-7/8,z-7/8")
    85      => SymOperation("x-3/4,y-1/4,z")
    86      => SymOperation("x-3/4,y-3/4,z-3/4")
    88      => SymOperation("x,y-3/4,z-7/8")
    125     => SymOperation("x-3/4,y-3/4,z")
    126     => SymOperation("x-3/4,y-3/4,z-3/4")
    129     => SymOperation("x-3/4,y-1/4,z")
    130     => SymOperation("x-3/4,y-1/4,z")
    133     => SymOperation("x-3/4,y-1/4,z-3/4")
    134     => SymOperation("x-3/4,y-1/4,z-3/4")
    137     => SymOperation("x-3/4,y-1/4,z-3/4")
    138     => SymOperation("x-3/4,y-1/4,z-3/4")
    141     => SymOperation("x,y-1/4,z-7/8")
    142     => SymOperation("x,y-1/4,z-7/8")
    146     => SymOperation("y,-x+y,z")
    148     => SymOperation("y,-x+y,z")
    151     => SymOperation("x,y,z-1/6")
    153     => SymOperation("x,y,z-5/6")
    154     => SymOperation("x,y,z-1/6")
    155     => SymOperation("y,-x+y,z")
    160     => SymOperation("y,-x+y,z")
    161     => SymOperation("y,-x+y,z")
    166     => SymOperation("y,-x+y,z")
    167     => SymOperation("y,-x+y,z")
    201     => SymOperation("x-3/4,y-3/4,z-3/4")
    203     => SymOperation("x-7/8,y-7/8,z-7/8")
    222     => SymOperation("x-3/4,y-3/4,z-3/4")
    224     => SymOperation("x-3/4,y-3/4,z-3/4")
])...)