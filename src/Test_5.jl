using OrdinaryDiffEq, Plots, DiffEqDevTools, LaTeXStrings

abstols = 1.0 ./ 10.0 .^ (3:7);  reltols = 1.0 ./ 10.0 .^ (3:7) 

M = zeros(3,3); M[1,1] = 1.0; M[2,2] = 1.0; u0 = [1.0, 1.0,-6.0]; tspan = [0.0, 1.0]

function f(dy, y, p, t)
  dy[1] = 0.5*y[2]^3*y[3]
  dy[2] = y[2]*y[3]/6.0
  dy[3] = y[3] + 6*y[1]/y[2]^3
end

function jac(J,u,p,t)
  J[:] = [0.0 1.5*1*(-6) 0.5;0.0 -1.0 1.0/6;6.0 -18.0 1]
  nothing
end

prob = ODEProblem(ODEFunction(f,mass_matrix = M), u0, tspan)

sol = solve(prob,Rodas5P(),abstol=1.0e-14,reltol=1.0e-14,maxiters=Int(1e8))
test_sol = TestSolution(sol)

prob = ODEProblem(ODEFunction(f,mass_matrix = M, jac = jac), u0, tspan)
setups = [Dict(:alg=>Velds4()),Dict(:alg=>ROS34PW2()),Dict(:alg=>ROS34PW3()),Dict(:alg=>ROS34PRw()),Dict(:alg=>Rodas23W()),Dict(:alg=>Rodas4P2()),
          Dict(:alg=>Rodas5P()),Dict(:alg=>Rodas5Pe()),Dict(:alg=>Rodas5Pr())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:l2,dense_errors=false,appxsol=test_sol,maxiters=Int(1e7),numruns=10)
plot(wp,xlabel=latexstring("\$l_2\$ error"))


