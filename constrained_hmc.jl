using HomotopyContinuation, DynamicPolynomials
using LinearAlgebra
using PyPlot
using DelimitedFiles

include("ellipse.jl")

# the stepsize \tau used in the proposal scheme
step_size = 0.2

# check the forward scheme and print information, if the flag equals 1.
check_rattle_flag = 1
# tolerance error
check_tol = 1e-6
backward_check_tol = 1e-6

# \xi: a k-dimensional vector
function xi(x)
  return [xi_i(x,i) for i in 1:k]
end

# \nabla\xi: a k\times d matrix
function grad_xi(x)
  tmp = zeros(k,d)
  for i in 1:k
    tmp[i,:] = grad_xi_i(x,i)
  end
  return tmp
end

# Hamiltonian energy 
function energy(x,v)
  return V(x) + norm(v) * 0.5
end

F_p = subs(F, p => p0)
# Compute all solutions for F_p, the starting system
result_p = solve(F_p)

# record the solutions
S_p0 = solutions(result_p)
num_sol_start_system = length(S_p0)
println("No. of real solutions for the starting systems: ", num_sol_start_system)

#Construct the PathTracker
tracker = pathtracker(F; parameters=p, generic_parameters=p0)

# solve algebraic equations for given parameters p
function find_solutions(p)
    # Create an empty array.
    S_p = similar(S_p0, 0)
    for s in S_p0
        result = track(tracker, s; target_parameters=p)
        # check that the tracking was successfull
       if is_success(result) && is_real(result)
         sol=solution(result)
	 push!(S_p, sol)
       end
    end
    return S_p
end

# compute several possible proposal states, at the current state (x,v)
function forward_rattle(x, v, step_size)
  grad_pot_vec = grad_V(x)
  # this should be a (k x d) matrix
  grad_xi_vec = grad_xi(x)
  coeff = - step_size * step_size * 0.5
  x_tmp = x + step_size * v + coeff * grad_pot_vec 
  p = vcat(x_tmp, step_size * grad_xi_vec[1,:])
  # find Lagrange multipliers for x
  S_p = find_solutions(p)
  n = length(S_p)
  # if we find at least one solutions
  if n > 0
    # randomly choose one
    j = rand(1:n)
    # compute the new state x^1
    x_1 = x_tmp + step_size * grad_xi_vec[1,:] * S_p[j][1].re
    if check_rattle_flag == 1 && norm(xi(x_1)) > check_tol
      println("x^1 is not on the level set, |xi(x^1)|=", norm(xi(x_1)))
      exit(1)
    end
    # prepare to compute the Lagrange multiplier lam_v
    v_tmp = v - 0.5 * step_size * (grad_pot_vec + grad_V(x_1)) + grad_xi_vec[1,:] * S_p[j][1].re
    grad_xi_vec_1 = grad_xi(x_1)
    mat_v_tmp = grad_xi_vec_1 * transpose(grad_xi_vec_1)
    # directly compute the Lagrange multiplier lam_v, by solving a linear system
    lam_v = - inv(mat_v_tmp) * grad_xi_vec_1 * v_tmp 
    # compute the updated velocity v^1
    v_1 = v_tmp + grad_xi_vec_1[1,:] * lam_v[1]
    return n, j, x_1, v_1
  end # no solutions are found, if we reach here
  return 0, 0, x, v
end

#backward check, similar to the forward_rattle function.
function backward_check(x1, v1, x, v, step_size)
  grad_pot_vec = grad_V(x1)
  # this should be a (k x d) matrix
  grad_xi_vec = grad_xi(x1)
  coeff = - step_size * step_size * 0.5
  x_tmp = x1 + step_size * v1 + coeff * grad_pot_vec 
  p = vcat(x_tmp, step_size * grad_xi_vec[1,:])
  # find Lagrange multipliers for x
  S_p = find_solutions(p)
  backward_found_flag = 0
  jj = 0
  # go through all solutions, and check if one solution is (x,v)
  for j in 1:length(S_p)
    # compute the new state x^2
    x_2 = x_tmp + step_size * grad_xi_vec[1,:] * S_p[j][1].re
    # first check whether the states are the same
    if norm(x_2 - x) > backward_check_tol
      continue
    else # if the states are the same, compute the velocity and check 
      # prepare to compute the Lagrange multiplier lam_v
      v_tmp = v1 - 0.5 * step_size * (grad_pot_vec + grad_V(x_2)) + grad_xi_vec[1,:] * S_p[j][1].re
      grad_xi_vec_2 = grad_xi(x_2)
      mat_v_tmp = grad_xi_vec_2 * transpose(grad_xi_vec_2)
      # directly compute the Lagrange multiplier lam_v, by solving a linear system
      lam_v = - inv(mat_v_tmp) * grad_xi_vec_2 * v_tmp 
      # compute the updated velocity v^2
      v_2 = v_tmp + grad_xi_vec_2[1,:] * lam_v[1]
      if norm(v_2 - v) < backward_check_tol # successful if we are here, both state and velocity are the same
        backward_found_flag = 1
	# record the index as well. It will be used to compute the M-H rate
	jj = j
	break 
      end
    end
  end
  return backward_found_flag, length(S_p), jj
end

function rand_draw_velocity(x)
  grad_xi_vec = grad_xi(x)
  # generate the orthnormal basis of the tangent space
  U_x = nullspace(grad_xi_vec)
  # generate normal Gaussian vector, as coefficients under the basis  
  coeff = randn(d-k)
  return U_x * coeff
end

# initial state
x0 = [3.0, 0.0]
v0 = [0.0, 1.0]

# total number of samples 
N = 50000

# array to store the samples 
sample_data = [zeros(2 * d) for i in 1:N]
forward_success_counter = 0
backward_success_counter = 0

# count the runtime 
@time begin
# the main loop 
for i in 1:N
  global x0, v0, forward_success_counter, backward_success_counter
  # first of all, randomly update the velocity 
  v0 = rand_draw_velocity(x0)
  # save the current state
  sample_data[i] = vcat(x0, v0)
  # compute proposal states
  n, j, x1, v1 = forward_rattle(x0, v0, step_size)
  if n > 0 # if one solution has been found, do backward check
    forward_success_counter += 1 
    # reverse the velocity, and do backward check
    found_flag, n_back, j_back = backward_check(x1, -v1, x0, -v0, step_size)
    if found_flag == 1 # if the backward check is passed 
      backward_success_counter += 1
      h = energy(x0, v0) 
      h_1 = energy(x1, v1) 
      # compute the MH-rate
      mh_rate = min(1, exp(h - h_1) * n / n_back)
      r = rand()
      if r < mh_rate # accept the proposal
        x0 = x1
        v0 = v1
      	continue
      end
    end
  end
# it is not necessary to flip the velocity, 
# because the velocity will be redrawn at the begining of the next step
#  v0 = -v0
end
# time end
end

println("counters: ", forward_success_counter, ' ', backward_success_counter)

# write the sampled states to the file
writedlm("./data/data.txt", sample_data)
