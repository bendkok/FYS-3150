include("ising.jl")

L = 20
s = ones(L,L)
sr = rand(-1:2:1, L,L)

J = 1
k = 1
T = 1
n = 2e5
#n = 0.25e5
#n = 0.05e5

Ising_model_test(J,s,n,10*T,k,"o")
#Ising_model_test(J,sr,n,T,k,"r")

#Ising_model_test(J,s,n,2.4,k,"o")
#Ising_model_test(J,sr,n,2.4,k,"r")



#E1,M1,s1,ye1,ym1 = Ising_model_test(J,s,n,3,k,"o")

"mE = mean_val(E)
mM = mean_val(M)
vE = varia(E)
vM = varia(M)

println(mE)
println(mM)
println(vE)
println(vM)"