using OrdinaryDiffEq, Plots, DiffEqDevTools, LaTeXStrings

setups = [Dict(:alg=>Rodas23W()),Dict(:alg=>Rodas3P()),Dict(:alg=>Rodas5Pe()),Dict(:alg=>Rodas5Pr())]

abstols = 1.0 ./ 10.0 .^ (2:10);  reltols = 1.0 ./ 10.0 .^ (2:10) 

function f!(du, u, fre, t)
    du[1] = u[1] - sin(fre*2*pi*t) + 1
    du[2] = 1
end

function f_analytic(u₀,p,t)
    [sin(fre*2*pi*t),t]
end

u0 = [0.0, 0.0]; tspan = (0.0, 1.0); 
M = zeros(2,2);  M[1,2] = 1; M[2,2] = 1; fre = 10.0
prob = ODEProblem(ODEFunction{true, SciMLBase.FullSpecialize}(f!, analytic = f_analytic, mass_matrix = M), u0, tspan, fre)

wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:L2, dense_errors=true,maxiters=Int(1e8),numruns=10)
plot(wp,xlabel=latexstring("\$L_2\$ error"))
