using OrdinaryDiffEq, Plots, DiffEqDevTools, LaTeXStrings, SparseArrays, Symbolics, LinearSolve

abstols = 1.0 ./ 10.0 .^ (2:10);  reltols = 1.0 ./ 10.0 .^ (2:10) 

nx = 500; dx = 2.0/nx; n = nx-1; x = zeros(n)
for i=1:n
  x[i] = -1.0 + i*dx
end
b1  = zeros(n); J = zeros(n,n)

b1[1]= -1.0/dx^2; b1[n]= 1.0/dx^2;
for i=1:n
  J[i,i] = -2/(dx^2);
  if i>1 J[i,i-1] = 1/(dx^2); end
  if i<n J[i,i+1] = 1/(dx^2); end
end
param = x, b1, J  

function f_analytic(u₀,p,t)
    x, b1, J = p  
    u = exp(t)*x.^3;
end

u0 = f_analytic(~,param,0.0); tspan = [0.0, 1.0]; 

function ode!(du,u,p,t)
    x, b1, J = p  
    g = exp(t)*(x.^3 -6*x - exp(t)*x.^6)
    rb = b1*exp(t)
    du[:] = J*u + u.^2 + rb + g;
end

du0 = copy(u0)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> ode!(du, u, param, 0.0), du0, u0)

setups = [Dict(:alg=>Scholz4_7(linsolve = KLUFactorization())),Dict(:alg=>ROS3P(linsolve = KLUFactorization())),Dict(:alg=>ROS3PR(linsolve = KLUFactorization())),Dict(:alg=>ROS3PRL(linsolve = KLUFactorization())),Dict(:alg=>ROS3PRL2(linsolve = KLUFactorization())),
          Dict(:alg=>ROS34PW1a(linsolve = KLUFactorization())),Dict(:alg=>ROS34PW1b(linsolve = KLUFactorization())),Dict(:alg=>ROS34PW2(linsolve = KLUFactorization())),Dict(:alg=>ROS34PW3(linsolve = KLUFactorization())),Dict(:alg=>ROS34PRw(linsolve = KLUFactorization())),
          Dict(:alg=>Rodas3P(linsolve = KLUFactorization())),Dict(:alg=>Rodas4P(linsolve = KLUFactorization())),Dict(:alg=>Rodas4P2(linsolve = KLUFactorization())),
          Dict(:alg=>Rodas5P(linsolve = KLUFactorization())),Dict(:alg=>Rodas5Pe(linsolve = KLUFactorization())),Dict(:alg=>Rodas5Pr(linsolve = KLUFactorization()))]
prob = ODEProblem(ODEFunction{true, SciMLBase.FullSpecialize}(ode!, analytic = f_analytic; jac_prototype = float.(jac_sparsity)), u0, tspan, param)
wp = WorkPrecisionSet(prob,abstols,reltols,setups; save_everystep=true, error_estimate=:l2,dense_errors=false,maxiters=Int(1e7),numruns=10)
plot(wp,xlabel=latexstring("\$l_2\$ error"))
