include("ising.jl")
using NPZ

function many_Ising_to_file(L,n)
    """
    A function that runs the Ising_model with many different values of T
    for a random LxL matrix s, and then writes the reulting E and M to a
    file.
    T goes between 2.0 and 2.6 with dT=0.05, ot not
    Takes the inputs L for length of matrix and n for number of cycles.
    """
    T = 2.0:0.05:2.6
    s = rand(-1:2:1, L,L)
    @loadbar "toFile" for t in T
        E0,M0,s0=Ising_model(1,s,n,t,1)
        E,M,s=Ising_model(1,s0,n/5,t,1)
        filewriter(E, "E_"*string(L)*"_"*string(t)*"_"*string(n)*".npy")
        filewriter(M, "M_"*string(L)*"_"*string(t)*"_"*string(n)*".npy")
    end
    println(L)
end

function filewriter(dataArray, filename = "test.npz")
    """
    A function that writes an array dataArray to a file filename.
    I had to include the full filepath be able to run the program in VSC.
    """
    path = "C:/Users/Bendik/Documents/FYS3150/Prosjekt 4/res/"
    f = path*filename
    npzwrite(f, dataArray)
    "vars = npzread(f)
    println(vars)"
end
#filewriter([1,2,2,3,4,5,3])

#many_Ising_to_file(20,10)
many_Ising_to_file(40,1e5)
many_Ising_to_file(60,1e5)
many_Ising_to_file(80,1e5)
many_Ising_to_file(100,1e5)


"""
[Running] julia "c:\\Users\\Bendik\\Documents\\FYS3150\\Prosjekt 4\\opp_e.jl"
40
60
80
100

[Done] exited with code=0 in 346.885 seconds


[Done] exited with code=0 in 2516.678 seconds

[Done] exited with code=0 in 2989.575 seconds

[Done] exited with code=0 in 7913.493 seconds

[Done] exited with code=0 in 2094.687 seconds

[Done] exited with code=0 in 1546.602 seconds

"""
