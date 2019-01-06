using Plots


# Plot of linked scalar sequences
x = 0:0.001:5
s1 = sin.(x)'
s2 = sin.(log.(range(1;length=5001, stop=exp(5))))'
l, m = linkseq(s1', s2', 1000, 50)

p = plot(x, [s1[:] s2[:]])
plot!(p, l, m; legend=false)
savefig(p, "docs/assets/2d_linked_seqs.svg")


# Plot of linked vector sequences
x = 0:0.001:5
s1 = [0.3; -0.7] * sin.(x)'
s2 = [0.7; 0.3] * sin.(log.(range(1;length=5001, stop=exp(5))))'
l, m = linkseq(s1, s2, 1000, 50)

p = plot(x, [s1[1,:] s2[1,:]], [s1[2,:] s2[2,:]])
plot!(p, l, m[1,:,:], m[2,:,:]; legend=false)
savefig(p, "docs/assets/3d_linked_seqs.svg")


# Plot different linked sequences
s1 = [0, 1, 2, 4, 5, 2, 3, 1, 0, 0]
s2 = [2, 4, 5, 6, 6, 5, 5, 4, 3, 1, 0, 0]
l, m = linkseq(s1, s2)
p = plot(s1)
plot!(p, s2)
plot!(p, l.+1, m; color=:black, legend=false)
savefig(p, "docs/assets/different_seqs.svg")
