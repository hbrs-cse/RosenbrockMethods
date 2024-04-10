using OrdinaryDiffEq, Plots

function f!(du, u, fre, t)
    du[1] = u[1] - sin(fre*2*pi*t) 
    du[2] = 0
end

u0 = [0.0, 0.0]; tspan = (0.0, 1.0); M = zeros(2,2);  M[2,2] = 1; fre = 10.0
fode = ODEFunction(f!, mass_matrix = M); 
prob = ODEProblem(fode, u0, tspan, fre)
sol = solve(prob, Rodas4(), dense = true)

P1 = scatter(sol.t,sol[1,:],label="Rodas4")
plot!(P1,sol,idxs=1,label="Interpolation")
tt = 0:0.001:tspan[2]
plot!(P1,tt,sin.(fre*2*pi*tt), label="solution",title = "Test 1", ylabel="y(t)")
