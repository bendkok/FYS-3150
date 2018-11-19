include("ising.jl")
L = 20
s = ones(L,L)

J = 1
k = 1
T0 = 1
T1 = 2.4
n = 30000 #seems fine

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
    plt.plot(Es,Pe)
    plt.grid()
    plt.xlabel("Energy")
    plt.ylabel("Probability")
    plt.savefig("C:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_d_"*string(n)*"_"*string(T)*"K.pdf")
    plt.show()
    return P,Pe
end

#P1,Pe1 = prob(s,J,k,T0,n,n)
P2,Pe2 = prob(s,J,k,T1,n*2,n)