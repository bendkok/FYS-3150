include("ising.jl")
import Printf

function do_stuff(n)
    s = [1 1; 1 1]
    J = 1
    k = 1
    #n = 1e7
    @time E,M,s0 = Ising_model(J,s,n,1,1)


    deg = [1,4,4,2,4,1] #the degeneracy
    Ei = [-8J, 0, 0, 8J, 0, -8J]
    Mi = [4, 2, 0, 0, -2, -4]
    T = 1
    Z, Em, Mm, CV, X = analytical_values(T, k, Ei, Mi, deg)

    Enum = mean_val(E)
    CVnum = varia(E)/(k*T^2)
    Mnum = mean_val(M)
    Xnum = varia(M)/(k*T)

    """println(Z)
    println()
    println(Em)
    println(Enum)
    println((1-Em/Enum)*1000)
    println()
    println(CV)
    println(CVnum)
    println((1-CV/CVnum)*1000)
    println()
    println(Mm)
    println(Mnum)
    println((1-Mm/Mnum)*1000)
    println()
    println(X)
    println(Xnum)
    println((1-X/Xnum)*1000)"""

    #s = @Printf.sprintf """\$Z\$ & %.2f  &	--- & --- \\\\ \\hline 
    #\$\\langle E \\rangle\$ & %.4f & %.4f & %.5f \\\\ \\hline
		#\$|M|\$ & %.4f & %.4f & %.5f \\\\ \\hline
		#\$C_V\$ & %.6f & %.6f & %.5f \\\\ \\hline
		#\$\\chi\$ & %.6f & %.6f & %.5f \\\\ \\hline
    #""", Z, Em, Enum, (1-Em/Enum)*1000, Mm, Mnum, (1-Mm/Mnum)*1000, CV, CVnum, (1-CV/CVnum)*1000, X, Xnum, (1-X/Xnum)*1000

    Printf.@printf("\$Z\$ & %.3f  &	--- & --- \\\\ \\hline\n",Z)
    Printf.@printf("\$\\langle E \\rangle\$ & %.5f & %.5f & %.6f \\\\ \\hline\n", Em, Enum, abs((1-Em/Enum)*1000))
    Printf.@printf("\$|M|\$ & %.6f & %.6f & %.6f \\\\ \\hline\n", Mm, Mnum, abs((1-Mm/Mnum)*1000))
    Printf.@printf("\$C_V\$ & %.6f & %.6f & %.6f \\\\ \\hline\n", CV, CVnum, abs((1-CV/CVnum)*1000))
    Printf.@printf("\$\\chi\$ & %.6f & %.6f & %.6f \\\\ \\hline\n", X, Xnum, abs((1-X/Xnum)*1000))
    println("")
end

println("n=1e4:")
do_stuff(1e4)
println("\nn=1e5:")
do_stuff(1e5)
println("\nn=1e6:")
do_stuff(1e6)
println("\nn=1e7:")
do_stuff(1e7)


"""


[Running] julia "c:\Users\Bendik\Documents\FYS3150\Prosjekt 4\opp_b.jl"
n=1e4:
  0.005056 seconds (40.03 k allocations: 4.240 MiB)
$Z$ & 5973.917  &	--- & --- \\ \hline
$\langle E \rangle$ & -7.98393 & -7.97920 & 0.592323 \\ \hline
$|M|$ & 3.994643 & 3.992201 & 0.611731 \\ \hline
$C_V$ & 0.128329 & 0.165967 & 226.779928 \\ \hline
$\chi$ & 0.016043 & 0.025939 & 381.516044 \\ \hline


n=1e5:
  0.078317 seconds (400.04 k allocations: 41.386 MiB, 4.67% gc time)
$Z$ & 5973.917  &	--- & --- \\ \hline
$\langle E \rangle$ & -7.98393 & -7.98320 & 0.091214 \\ \hline
$|M|$ & 3.994643 & 3.994360 & 0.070819 \\ \hline
$C_V$ & 0.128329 & 0.134118 & 43.159349 \\ \hline
$\chi$ & 0.016043 & 0.017008 & 56.751048 \\ \hline


n=1e6:
  0.608158 seconds (4.00 M allocations: 377.842 MiB, 12.72% gc time)
$Z$ & 5973.917  &	--- & --- \\ \hline
$\langle E \rangle$ & -7.98393 & -7.98398 & 0.006973 \\ \hline
$|M|$ & 3.994643 & 3.994664 & 0.005276 \\ \hline
$C_V$ & 0.128329 & 0.127871 & 3.580466 \\ \hline
$\chi$ & 0.016043 & 0.015972 & 4.472392 \\ \hline


n=1e7:
  5.761999 seconds (40.00 M allocations: 3.655 GiB, 7.78% gc time)
$Z$ & 5973.917  &	--- & --- \\ \hline
$\langle E \rangle$ & -7.98393 & -7.98383 & 0.012769 \\ \hline
$|M|$ & 3.994643 & 3.994612 & 0.007643 \\ \hline
$C_V$ & 0.128329 & 0.129153 & 6.376069 \\ \hline
$\chi$ & 0.016043 & 0.016125 & 5.061573 \\ \hline


[Done] exited with code=0 in 14.469 seconds



"""
