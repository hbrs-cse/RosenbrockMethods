using OrdinaryDiffEq, Plots, StaticArrays

methods = Rodas23W(), Rodas3P(), Rodas5Pe(), Rodas5Pr(), Rodas5P()

function f!(du, u, fre, t) #-- inplace problem
    du[1] = u[1] - sin(fre*2*pi*t) 
    du[2] = 0
end

function f(u, fre, t) #-- out of place problem
    du = [u[1] - sin(fre*2*pi*t), 0]
end

function fsa(u, fre, t) #-- static array problem
    SA[u[1] - sin(fre*2*pi*t), 0]
end

function f22!(du, u, fre, t) #-- matrix problem
    du[1,1] = u[1,1] - sin(fre*2*pi*t) 
    du[2,1] = 0
    du[1,2] = 0
    du[2,2] = 0
end

function fs(u, fre, t) #-- scalar problem
    du = fre*2*pi*cos(fre*2*pi*t)
end

condition(u, t, integrator) = t - 1.0 
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

function test1(method) #-- inplace + stiff_add_step
  u0 = [0.0, 0.0]; tspan = (0.0, 2.0); M = zeros(2,2);  M[2,2] = 1; fre = 10.0
  fode = ODEFunction(f!, mass_matrix = M); prob = ODEProblem(fode, u0, tspan, fre) 
  sol = solve(prob, method, dense = true, callback = cb)
  tt = LinRange(0.0,1.0,100)
  err = maximum(abs.(sol(tt; idxs = 1) - sin.(fre*2*pi*tt)))
  println(method," ",sol.destats," ",err)
end

function test2(method) #-- out of place
  u0 = [0.0, 0.0]; tspan = (0.0, 1.0); M = zeros(2,2);  M[2,2] = 1; fre = 10.0
  fode = ODEFunction(f, mass_matrix = M); prob = ODEProblem(fode, u0, tspan, fre) 
  sol = solve(prob, method, dense = true)
  tt = LinRange(0.0,1.0,100)
  err = maximum(abs.(sol(tt; idxs = 1) - sin.(fre*2*pi*tt)))
  println(method," ",sol.destats," ",err)
end
   
function test3(method) #-- static array
  u0s = SA[0.0, 0.0];  tspan = (0.0, 1.0); M = zeros(2,2);  M[2,2] = 1; fre = 10.0
  fode = ODEFunction(fsa, mass_matrix = M); prob = ODEProblem(fode, u0s, tspan, fre) 
  sol = solve(prob, method, dense = true)
  tt = LinRange(0.0,1.0,100)
  err = maximum(abs.(sol(tt; idxs = 1) - sin.(fre*2*pi*tt)))
  println(method," ",sol.destats," ",err)
end

function test4(method) #-- matrix problem
  u0 = [0.0 0.0;0.0 0.0]; tspan = (0.0, 1.0); M = zeros(4,4);  M[2,2] = 1;M[3,3] = 1;M[4,4] = 1; fre = 10.0
  fode = ODEFunction(f22!, mass_matrix = M); prob = ODEProblem(fode, u0, tspan, fre)
  sol = solve(prob, method, dense = true)
  tt = LinRange(0.0,1.0,100)
  err = maximum(abs.(sol(tt; idxs = 1) - sin.(fre*2*pi*tt)))
  println(method," ",sol.destats," ",err)
end

function test5(method) #-- scalar problem + stiff_add_step
  u0 = 0.0; tspan = (0.0, 1.0); fre = 10.0
  fode = ODEFunction(fs); prob = ODEProblem(fode, u0, tspan, fre)
  sol = solve(prob, method, dense = true, callback = cb)
  tt = LinRange(0.0,1.0,100)
  err = maximum(abs.(sol(tt; idxs = 1) - sin.(fre*2*pi*tt)))
  println(method," ",sol.destats," ",err)
end

for method in methods test1(method); end
for method in methods test2(method); end
for method in methods test3(method); end
for method in methods test4(method); end
for method in methods test5(method); end
