include("ising.jl")

s = [[1 1]; [1 1]]
J = 1
k = 1
E,M,s0 = Ising_model(J,s,40000,1,1)


deg = [1,4,4,2,4,1] #the degeneracy
Ei = [-8J, 0, 0, 8J, 0, -8J]
Mi = [4, 2, 0, 0, -2, -4]
T = 1
Z, Em, Mm, CV, X = analytical_values(T, k, Ei, Mi, deg)


println(Z)
println()
println(Em)
println(mean_val(E))
println()
println(CV)
println(varia(E)/(k*T^2))
println()
println(Mm)
println(mean_val(M))
println()
println(X)
println(varia(M)/(k*T))
