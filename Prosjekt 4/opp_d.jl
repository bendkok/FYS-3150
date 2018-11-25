include("ising.jl")
L = 20
#s = ones(L,L)
s = rand(-1:2:1, L,L)

J = 1
k = 1
T0 = 1
T1 = 2.4
n = 1e5 #seems fine

function prob(s,J,k,T,n,n1)
    """
    A function that finds the probability of each E
    """
    #first we reach the steady state situation
    E0,M0,s0 = Ising_model(J,s,n,T,k)
    #then compute for reals 
    E,M,s = Ising_model(J,s0,n1,T,k)
    #array with the amounts
    P = zeros(Int64, Int64(maximum(E)-minimum(E)+1)) 
    for k in E #goes 
        P[Int64(k-minimum(E))+1] += 1
    end
    Es = minimum(E):maximum(E) #array with the possible energy-levels
    Pe = P/length(E) #array with the probabilities
    println(sum(Pe))
    println(sum(P))
    println(minimum(E))
    println(maximum(E))
    println(maximum(E)-minimum(E))
    println(varia(E))
    println(sqrt(varia(E)))
    println()
    plt.plot(Es,Pe)
    plt.grid()
    plt.xlabel("Energy")
    plt.ylabel("Probability")
    plt.savefig("C:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_d_"*string(Int64(n1))*"_"*string(T*10)*"K.pdf")
    plt.show()
    return P,Pe
end

P1,Pe1 = prob(s,J,k,T0,n,n)
P2,Pe2 = prob(s,J,k,T1,n,2*n)

#=

[Running] julia "c:\Users\Bendik\Documents\FYS3150\Prosjekt 4\opp_d.jl"
1.0
20001
-800
-772
28
9.435921963901809
3.071794583610989

1.0
10001
-680
-272
408
3177.9853509049094
56.373622829342


[Done] exited with code=0 in 31.652 seconds


[Running] julia "c:\Users\Bendik\Documents\FYS3150\Prosjekt 4\opp_d.jl"
1.0
100001
-800
-772
28
9.097688543114575
3.0162374812197026

1.0
200001
-716
-288
428
3241.8300982407086
56.93707138798683


[Done] exited with code=0 in 99.655 seconds

=#
