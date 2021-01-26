# --- SCALARS ---
const DEFAULT_ATOL = 1e-12        # absolute tolerance for approximate equality
const NULL_ATOL    = 1e-11        # absolute tolerance for nullspace
const MAX_SGNUM    = (2, 17, 230) # number of space groups in dimensions 1, 2, and 3

# --- VECTORS ---
# arbitrary test vector for e.g. evaluating KVecs lines/planes/volumes; 3D by default
const TEST_αβγ = [0.123123123, 0.456456456, 0.789789789]

# --- LISTS ---
# The following space group numbers are symmorphic; in 3D they each embody one of 
# the 73 arithmetic crystal classes. Basically, this is the combination of the 14  
# Bravais lattices with the 32 point groups/geometric crystal classes. Obtained from e.g.
#   > tuple((tuple([i for i in 1:MAX_SGNUM[D] if issymmorph(i, D)]...) for D = 1:3)...)
const SYMMORPH_SGNUMS = (# 1D (1st index; all elements are symmorph)
                         (1, 2), 
                         # 2D (2nd index)
                         ((1, 2, 3, 5, 6, 9, 10, 11, 13, 14, 15, 16, 17)), 
                         # 3D (3rd index)
                         (1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47,65,
                          69,71,75,79,81,82,83,87,89,97,99,107,111,115,119,121,
                          123,139,143,146,147,148,149,150,155,156,157,160,162,
                          164,166,168,174,175,177,183,187,189,191,195,196,197,
                          200,202,204,207,209,211,215,216,217,221,225,229)
                        )