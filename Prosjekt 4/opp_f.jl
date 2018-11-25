include("ising.jl")
using NPZ

k = 1

function read_from_file(L,t,n)
    f = "C:/Users/Bendik/Documents/FYS3150/Prosjekt 4/res/E_"*string(L)*"_"*string(t)*"_"*string(n)*".npy"
    m = "C:/Users/Bendik/Documents/FYS3150/Prosjekt 4/res/M_"*string(L)*"_"*string(t)*"_"*string(n)*".npy"
    vars0 = npzread(f)
    vars1 = npzread(m)
    return vars0,vars1
end

function find_values(L,n,T)
    Em = []
    Mm = []
    Cv = []
    Xi = []

    for t in T
        e,m = read_from_file(L,t,n)
        
        em = mean_val(e)
        mm = mean_val(m)
        cv = varia(e)/(k*t^2)
        xi = varia(m)/(k*t)
        
        push!(Em, em)
        push!(Mm, mm)
        push!(Cv, cv)
        push!(Xi, xi)
    end
    return Em,Mm,Cv,Xi
end


M = []
Xi = []
T = 2.0:0.01:2.3
for t in T
    E40, M40, CV40, Xi40 = find_values(40,1e5,t)
    #E60, M60, CV60, Xi60 = find_values(40,1e5,t)
    #E80, M80, CV80, Xi80 = find_values(40,1e5,t)
    #E100, M100, CV100, Xi100 = find_values(40,1e5,t)

    push!(M, mean_val(M40))
    push!(Xi, mean_val(Xi40))
end
Mt = zeros(length(T))
Xit = zeros(length(T))
for i = 1:length(T)
    Mt[i] = abs(T[i] - M[i]^8)
    Xit[i] = abs(T[i] - Xi[i]^(4/7))
end

plt.plot(T,Mt)
plt.grid()
plt.show()
plt.plot(T,Xit)
plt.grid()
plt.show()

