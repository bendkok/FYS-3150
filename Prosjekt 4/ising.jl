import PyPlot
const plt = PyPlot
using Statistics

function energy(J, s)
    """
    A function that finds the energy for one cycle of the Monte Carlo
    method/the Ising model.
    
    J is a coupling constant expressing the strength of the interaction
    between neighboring spins, and s is a 2D-matrix with spins.

    Outputs E, the sum of the energies.
    """
    N = length(s)
    L = Int64(sqrt(N))
    E = 0 #the total energy
    even = L%2==0 #if L is an even number
    #goes through the all the spins
    for l = 1:L
        if l%2==0 #we have to make sure to not double count the pairs
            p = 2 #so when l is even we take the even values of k,
        else      #and when it is odd we take the odd values
            p = 1
        end
        for k = p:2:L
            if l == 1 #if the value is in the left column
                E += s[l,k]*s[L,k]
            else
                E += s[l,k]*s[l-1,k]
            end
            if l == L #if the value is in the rigth column
                if even #if L is odd we have allready counted it
                    E += s[l,k]*s[1,k]
                end
            else
                E += s[l,k]*s[l+1,k]
            end
            if k == 1 #if the value is in the top column
                E += s[l,k]*s[l,L]
            else
                E += s[l,k]*s[l,k-1]
            end
            if k == L #if the value is in the bottom column
                if even #if L is odd we have allready counted it
                    E += s[l,k]*s[l,1]
                end
            else
                E += s[l,k]*s[l,k+1]
            end
        end
    end
    return -J*E
end

function magnetization(s)
    """
    A function that finds the magnetization for a matrix of spins s.
    Outputs M, the sum of the spins.
    """
    M = sum(s)
    return M
end

function mean_val(x)
    """
    A function that finds the mean value for an array x
    """
    return Statistics.mean(x)
end

function varia(x)
    """
    A function that finds the variance for an array x
    """
    return Statistics.var(x)
end

function delta_E(J, s, l, k)
    """
    A function that finds ΔE for s[l,k]
    """
    dE = 0
    L = Int64(sqrt(length(s)))
    if l == 1 #if the value is in the left column
        dE += s[L,k]
    else
        dE += s[l-1,k]
    end
    if l == L #if the value is in the rigth column
        dE += s[1,k]
    else
        dE += s[l+1,k]
    end
    if k == 1 #if the value is in the top column
        dE += s[l,L]
    else
        dE += s[l,k-1]
    end
    if k == L #if the value is in the bottom column
        dE += s[l,1]
    else
        dE += s[l,k+1]
    end
    return dE*2*J*s[l,k]
end

function Ising_model(J,s,n,T,k)
    """
    A function that solves the Ising model for a matrix of spins s.
    Takes the inputs for the constants J and k, as well as the 
    temperature T and the number of cycles s.
    """
    N = length(s)
    L = Int64(sqrt(N))
    E = [0] #a list to hold all the energies
    M = [0] #a list to hold all the magnetizations
    beta = 1/(T*k)
    
    e1 = energy(J,s) #the values for the initial spins
    m1 = magnetization(s)
    "push!(E, e1)
    push!(M, m1)"
    E[1] = e1
    M[1] = m1

    rand_koor = zeros(L,2)
    for i = 1:n 
        #finn L^2 random koordinater
        for l = 1:L
            rand_koor[l,:] = rand(1:L, 2)
        end
        #gå gjennom alle koordinater
        for j = 1:length(L) #rand_koor[:,1]
            l = Int64(rand_koor[j,1])
            k = Int64(rand_koor[j,2])

            #flip spin
            s[l,k]*=-1
            #find dE
            dE = delta_E(J, s, l, k)
            #test dE <= 0
            if dE > 0 #false: keep fliped
                #få random tall r e[0,1]
                r = rand()
                w = exp(-dE*beta)
                #test r <= w
                if r[1] > w #false: keep fliped
                    s[l,k]*=-1
                end
            end
        end
        #gjør resten:

        #finds the energy and magnetization for current spins 
        e = energy(J,s) 
        m = magnetization(s)
        push!(E, Int64(e))
        push!(M, Int64(m))
    end
    return E,M,s
end

function Ising_model_test(J,s,n_max,T,k,r_o)
    """
    A function that test how many Monte Carlo cycles that are needed
    to reach equlibrium by plotting the mean E and M against number of
    cycles.
    """
    E,M,s0 = Ising_model(J,s,n_max,T,k) #finds E and M
    ye = zeros(length(E))
    ym = zeros(length(E))
    for i = 1:length(E) #finds the mean E and M for the current value
        ye[i] = mean_val(E[1:i]) #and all the one preciding
        ym[i] = abs(mean_val(M[1:i]))
    end
    x = 1:length(E)
    #plots the results
    plt.plot(x,ye)
    plt.grid()
    plt.xlabel("Cycles")
    plt.ylabel("Mean energy")
    plt.savefig("C:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_c_e_"*string(T)*"K_"*string(n_max)*"_"*r_o*".pdf")
    plt.show()

    plt.plot(x,ym)
    plt.grid()
    plt.xlabel("Cycles")
    plt.ylabel("Mean magnetization")
    plt.savefig("C:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_c_m_"*string(T)*"K_"*string(n_max)*"_"*r_o*".pdf")
    plt.show()
    
    return E,M,s,ye,ym
end

function analytical_values(T, k, Ei, Mi, deg)
    """
    A function that finds the analytical values for the partition 
    function Z, expectation values for the energy Em, the mean 
    absolute value of the magnetic moment Mm, the specific heat CV,
    and the susceptibility X

    Takes the inputs temperature T, 
    """
    beta = 1/(T*k)
    Z = 0 #the partition function
    for i = 1:length(Ei) 
       Z += deg[i]*exp(-beta*Ei[i]) 
    end

    Em = 0 #expectation values for the energy
    for i = 1:length(Ei) 
        Em += deg[i]*Ei[i]*exp(-beta*Ei[i])
    end
    Em *= 1/Z

    Mm = 0 #the mean absolute value of the magnetic moment
    for i = 1:length(Mi) 
        Mm += deg[i]*abs(Mi[i])*exp(-beta*Ei[i])
    end
    Mm = abs(Mm/Z)

    E2 = 0 #expectation values for the energy^2
    for i = 1:length(Ei) 
        E2 += deg[i]*Ei[i]^2*exp(-beta*Ei[i])
    end
    E2 *= 1/Z
    CV = (E2 - Em^2)/(k*T^2) #the specific heat 

    M2 = 0 #the mean absolute value of the magnetic moment^2
    for i = 1:length(Mi) 
        M2 += deg[i]*Mi[i]^2*exp(-beta*Ei[i])
    end
    M2 *= 1/Z
    X = (M2 - Mm^2)*beta #the susceptibility 
    
    return Z, Em, Mm, CV, X
end
