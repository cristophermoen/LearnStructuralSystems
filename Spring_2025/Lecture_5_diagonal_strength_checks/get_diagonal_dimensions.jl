t = 5/64
D = 1 + 11/64
B = 1 + 23/32


cross_section = [(t/2, D), (t/2, t/2), (B-t/2, t/2), (B-t/2, D)]

x = [cross_section[i][1] for i in eachindex(cross_section)]
y = [cross_section[i][2] for i in eachindex(cross_section)]