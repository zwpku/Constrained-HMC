using HomotopyContinuation 
using LinearAlgebra

d = 3
k = 1

R = 1.0
r = 0.5 

@polyvar lam[1:k] p[1:6]  
F = [(R^2 - r^2 + (p[1] + p[4] * lam[1])^2 + (p[2] + p[5] * lam[1])^2 + (p[3] + p[6] * lam[1])^2)^2 - 4 * R^2 * ((p[1] + p[4] * lam[1])^2 + (p[2] + p[5] * lam[1])^2)]						
#Randomly choose a start system
p0 = randn(Complex{Float32}, 6)
F0 = subs(F, p => p0)
res0 = solve(F0)
S0 = solutions(res0)

#Construct the tracker
tracker = pathtracker(F; parameters=p, generic_parameters=p0)

#Let's solve another system using PathTracking
p_current = randn(6)
result = track(tracker, S0[1]; target_parameters=p_current, accuracy=1e-9)
sol=solution(result)

println("Result=", result)

println("Accuracy=", accuracy(result), "\tResidual=", residual(result))

println("However, F(sol)=", F[1](lam=>sol, p=>p_current))
