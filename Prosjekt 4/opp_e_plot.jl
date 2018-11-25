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

function plot_helper(x,T,leg,ylab,filename)
    for i in x
        plt.plot(T[2:length(T)],i[2:length(T)]/1) 
    end
    plt.grid()
    plt.xlabel("Temperature T")
    plt.ylabel(ylab*"L^2")
    plt.legend(leg, loc="best")
    plt.savefig("C:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_e_"*filename*".pdf")
    plt.show()
end

function plot_stuff(Ls,T,n)
    e = []
    m = []
    cv = []
    xi = []
    for L in Ls
        Em,Mm,Cv,Xi = find_values(L,n,T)
        push!(e,Em/L^2)
        push!(m,Mm/L^2)
        push!(cv,Cv/L^2)
        push!(xi,Xi/L^2)
    end
    plot_helper(e,T,Ls,"Mean energy <E>/","me0")
    plot_helper(m,T,Ls,"Mean magnetization |M|/","mm0")
    plot_helper(cv,T,Ls,"Specific heat C_V/","cv0")
    plot_helper(xi,T,Ls,"Susceptibility X/","xi0")
end

#read_from_file(40,2.0,60000)
#println(find_values(40,60000))
plot_stuff([40,60,80,100],2.0:0.05:2.6,1e5)
#plot_stuff([40,60,80,100],2.0:0.05:2.3,1e5)
