using Crystalline

# print reality induced degeneracy properties of all irreps
linelen = 60
for sgnum = 1:230
    println("\n\nSpace group $(sgnum)\n", "──┬", '─'^linelen)
    for kidx = 1:length(LGIRS[sgnum])
        realify(LGIRS[sgnum][kidx], true)
        if kidx ≠ length(LGIRS[sgnum]) 
            println("──┼", '─'^linelen)
        end
    end
    println("──┴", '─'^linelen)
end