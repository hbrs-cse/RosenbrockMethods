using OrdinaryDiffEq, Plots, DiffEqDevTools, LaTeXStrings

abstols = 1.0 ./ 10.0 .^ (2:10);  reltols = 1.0 ./ 10.0 .^ (2:10) 

function f!(du, u, p, t)
    du[1] = -u[1]
    du[2] = u[2] - (1 - t^2)^4
end

function f_analytic(u₀,p,t)
    [exp(-t),(1-t^2)^4]
end
 
tspan = [0.0, 10.0]; u0 = [1.0,1.0]
M = zeros(2,2);  M[1,1] = 1;
prob = ODEProblem(ODEFunction{true, SciMLBase.FullSpecialize}(f!, analytic = f_analytic, mass_matrix = M), u0, tspan)

setups = [Dict(:alg=>Scholz4_7()),Dict(:alg=>ROS2S()),Dict(:alg=>ROS2PR()),Dict(:alg=>ROS3()),Dict(:alg=>ROS3P()),Dict(:alg=>ROS3PR()),Dict(:alg=>ROS3PRL()),Dict(:alg=>ROS3PRL2()),
          Dict(:alg=>ROS34PW1a()),Dict(:alg=>ROS34PW1b()),Dict(:alg=>ROS34PW2()),Dict(:alg=>ROS34PW3()),Dict(:alg=>ROS34PRw()),Dict(:alg=>Rodas5P())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:L2, dense_errors=true,maxiters=Int(1e8),numruns=10)
P1 = plot(wp,xlabel=latexstring("\$L_2\$ error"))

setups = [Dict(:alg=>ROS34PRw()),Dict(:alg=>Rodas23W()),Dict(:alg=>Rodas3()),Dict(:alg=>Rodas3P()),Dict(:alg=>Rodas4()),Dict(:alg=>Rodas4P()),
          Dict(:alg=>Rodas5P()),Dict(:alg=>Rodas5Pe()),Dict(:alg=>Rodas5Pr())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:L2, dense_errors=true,maxiters=Int(1e8),numruns=10)
P2 = plot(wp,xlabel=latexstring("\$L_2\$ error"))

plot(P1,P2)