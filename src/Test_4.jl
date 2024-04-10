using OrdinaryDiffEq, Plots, DiffEqDevTools, LaTeXStrings
gr()
default(titlefont = 6, legendfontsize = 4, guidefont = 5, tickfont = 5, markersize = 3.5)

abstols = 1.0 ./ 10.0 .^ (2:10);  reltols = 1.0 ./ 10.0 .^ (2:10) 

tspan = [0.0, 10.0]; u0 = [0.0, 1.0, 0.0]
J = [-0.8 12.5 0; -12.5 -0.8 0; 0 0 -100]

function f(du,u,J,t)
  du .= J*u 
end
 
prob = ODEProblem(ODEFunction(f), u0, tspan, J)

sol = solve(prob,Rodas5P(), abstol = 1.0e-14, reltol = 1.0e-14, maxiters = Int(1e8))
test_sol = TestSolution(sol)

setups = [Dict(:alg=>GRK4T()),Dict(:alg=>GRK4A()),Dict(:alg=>Veldd4()),Dict(:alg=>Veldd4()),Dict(:alg=>Velds4()),Dict(:alg=>Scholz4_7()),Dict(:alg=>Ros4LStab()),
         Dict(:alg=>Rodas4()),Dict(:alg=>Rodas42()),Dict(:alg=>Rosenbrock23()),Dict(:alg=>Rosenbrock32())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:l2,dense_errors=false,appxsol=test_sol,maxiters=Int(1e8),numruns=10)
P1 = plot(wp,xlabel=latexstring("\$l_2\$ error"))

setups = [Dict(:alg=>ROS2S()),Dict(:alg=>ROS2PR()),Dict(:alg=>ROS3()),Dict(:alg=>ROS3P()),Dict(:alg=>ROS3PR()),Dict(:alg=>ROS3PRL()),Dict(:alg=>ROS3PRL2()),
          Dict(:alg=>ROS34PW1a()),Dict(:alg=>ROS34PW1b()),Dict(:alg=>ROS34PW2()),Dict(:alg=>ROS34PW3()),Dict(:alg=>ROS34PRw()),Dict(:alg=>Rodas4())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:l2,dense_errors=false,appxsol=test_sol,maxiters=Int(1e8),numruns=10)
P2 = plot(wp,xlabel=latexstring("\$l_2\$ error"))

setups = [Dict(:alg=>Rosenbrock23()),Dict(:alg=>Rodas23W()),Dict(:alg=>Rodas3()),Dict(:alg=>Rodas3P()),Dict(:alg=>Rodas4()),Dict(:alg=>Rodas4P()),
          Dict(:alg=>Rodas4P2()),Dict(:alg=>Rodas5P()),Dict(:alg=>Rodas5Pe()),Dict(:alg=>Rodas5Pr())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:l2,dense_errors=false,appxsol=test_sol,maxiters=Int(1e8),numruns=10)
P3 = plot(wp,xlabel=latexstring("\$l_2\$ error"))

plot(P1,P2,P3)


