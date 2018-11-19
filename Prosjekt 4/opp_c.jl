include("ising.jl")

L = 20
s = ones(L,L)

J = 1
k = 1
T = 1
n = 40000

#E0,M0,s0,ye0,ym0 = Ising_model_test(J,s,n,T,k,"o")
E1,M1,s1,ye1,ym1 = Ising_model_test(J,s,n,2.4,k,"o")

sr = rand(-1:2:1, L,L)
#E2,M2,s2,ye2,ym2 = Ising_model_test(J,sr,n,T,k,"r")
E3,M3,s3,ye3,ym3 = Ising_model_test(J,sr,n,2.4,k,"r")

"mE = mean_val(E)
mM = mean_val(M)
vE = varia(E)
vM = varia(M)

println(mE)
println(mM)
println(vE)
println(vM)"