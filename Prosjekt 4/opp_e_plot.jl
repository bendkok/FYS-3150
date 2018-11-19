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

function plot_stuff(Ls,T)
    for L in Ls
        
    end
    """plt.plot(x,ye)
    plt.grid()
    plt.xlabel("Temperature T")
    plt.ylabel("Mean energy")
    plt.savefig("C:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_e.pdf")
    plt.show()"""
end

#read_from_file(40,2.0,60000)
println(find_values(40,60000))
