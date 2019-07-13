using HomotopyContinuation, DynamicPolynomials
using LinearAlgebra
using PyPlot
using DelimitedFiles
using Printf

include("ellipse.jl")

# the stepsize \tau used in the proposal scheme
step_size = 0.5

# check the forward scheme and print information, if the flag equals 1.
check_rattle_flag = 1
# tolerance error
check_tol = 1e-6
backward_check_tol = 1e-6
new_sol_tol = 1e-6

# the code will be slower, without PathTracking
path_tracking_flag = 1

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
  return V(x) + dot(v,v) * 0.5
end

# prepare the start system
if path_tracking_flag == 1
  F_p = subs(F, p => p0)
  # Compute all solutions for F_p, the starting system
  # according to the package's usage, Total Degree Homotopy is used.
  result_p = solve(F_p)

  # record the solutions
  S_p0 = solutions(result_p)
  num_sol_start_system = length(S_p0)

  @printf("Starting systems: no. of real solutions = %d\n", num_sol_start_system)

  #Construct the PathTracker
  tracker = pathtracker(F; parameters=p, generic_parameters=p0)
end

# solve algebraic equations for given parameters p, with path tracking
function find_solutions_by_tracking(p)
    # Create an empty array.
    S_p = similar(S_p0, 0)
    for s in S_p0
        result = track(tracker, s; target_parameters=p)
        # check that the tracking was successfull
       if is_success(result) && is_real(result)
         sol=solution(result)
	 new_sol_flag = 1
	 for i in 1:length(S_p)
	    if euclidean_distance(S_p[i], sol) < new_sol_tol
	      new_sol_flag = 0
	      break
	    end
	 end
	 if new_sol_flag == 1
	   push!(S_p, sol)
	 end
       end
    end
    return S_p
end

# solve equations without path tracking
function find_solutions_total_degreee(p_current)
  F_p = subs(F, p => p_current)
  # Compute all solutions for F_p  
  # according to the package's usage, Total Degree Homotopy is used.
  result_p = solve(F_p)
  # record the solutions
  S_p = solutions(result_p; only_real=true)
  return S_p
end

function find_solutions(p)
  if path_tracking_flag == 1
    S_p = find_solutions_by_tracking(p)
  else 
    S_p = find_solutions_total_degreee(p)
  end
  n = length(S_p)
  lambda_vec = zeros(k,n)
  # extract the real part
  for i in 1:n
    lambda_vec[:,i] = [S_p[i][j].re for j in 1:k]
  end
  return lambda_vec
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
  lam_x = find_solutions(p)
  n = length(lam_x)
  # if we find at least one solutions
  if n > 0
    dist = zeros(n)
    # compute the distance wrt x for each proposal
    for i in 1:n
      dist[i] = norm(x_tmp + step_size * transpose(grad_xi_vec) * lam_x[:,i] - x)
    end
    # sort accrording to the distance
    perm = sortperm(dist)
    # randomly choose one
    j = rand(1:n)
    # compute the new state x^1
    x_1 = x_tmp + step_size * transpose(grad_xi_vec) * lam_x[:,j]
    if check_rattle_flag == 1 && norm(xi(x_1)) > check_tol
      println("x^1 is not on the level set, |xi(x^1)|=", norm(xi(x_1)))
      exit(1)
    end
    # prepare to compute the Lagrange multiplier lam_v
    v_tmp = v - 0.5 * step_size * (grad_pot_vec + grad_V(x_1)) + transpose(grad_xi_vec) * lam_x[:,j] 
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
  lam_x = find_solutions(p)
  backward_found_flag = 0
  jj = 0
  # go through all solutions, and check if one solution is (x,v)
  for j in 1:length(lam_x)
    # compute the new state x^2
    x_2 = x_tmp + step_size * transpose(grad_xi_vec) * lam_x[:,j]
    # first check whether the states are the same
    if norm(x_2 - x) > backward_check_tol
      continue
    else # if the states are the same, compute the velocity and check 
      # prepare to compute the Lagrange multiplier lam_v
      v_tmp = v1 - 0.5 * step_size * (grad_pot_vec + grad_V(x_2)) + transpose(grad_xi_vec) * lam_x[:,j] 
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
      else 
        printfln("the same state, different velocity!")
	println(x1, v1, x, v)
      end
    end
  end
  return backward_found_flag, length(lam_x), jj
end

function rand_draw_velocity(x)
  grad_xi_vec = grad_xi(x)
  # generate the orthnormal basis of the tangent space
  U_x = nullspace(grad_xi_vec)
  # generate normal Gaussian vector, as coefficients under the basis  
  coeff = randn(d-k)
  return U_x * coeff
end

# total number of samples 
N = 50000

# array to store the samples 
sample_data = [zeros(2 * d) for i in 1:N]
forward_success_counter = 0
backward_success_counter = 0
stat_success_counter = 0

# upper bound of solution number, let us assume there are at most 10 solutions
max_no_sol = 10
stat_num_of_solution_forward = zeros(max_no_sol+1)
stat_num_of_solution_backward = zeros(max_no_sol+1)

# count the runtime 
@time begin
# the main loop 
for i in 1:N
  global x0, v0, forward_success_counter, backward_success_counter, stat_success_counter, stat_num_of_solution_forward, stat_num_of_solution_backward
  # first of all, randomly update the velocity 
  v0 = rand_draw_velocity(x0)
  # save the current state
  sample_data[i] = vcat(x0, v0)
  # compute proposal states
  n, j, x1, v1 = forward_rattle(x0, v0, step_size)
  if n <= max_no_sol 
    stat_num_of_solution_forward[n+1] += 1
  else 
    @printf("Warning: No. of solutions in forward rattle (=%d) is larger than upper bound (=%d)!", n, max_no_sol)
  end
  if n > 0 # if one solution has been found, do backward check
    forward_success_counter += 1 
    # reverse the velocity, and do backward check
    found_flag, n_back, j_back = backward_check(x1, -v1, x0, -v0, step_size)
    if n_back <= max_no_sol 
      stat_num_of_solution_backward[n_back+1] += 1
    else 
      @printf("Warning: No. of solutions in backward check (=%d) is larger than upper bound (=%d)!", n_back, max_no_sol)
    end
    if found_flag == 1 # if the backward check is passed 
      backward_success_counter += 1
      h = energy(x0, v0) 
      h_1 = energy(x1, v1) 
      # compute the MH-rate
      mh_rate = min(1, exp(h - h_1) * n / n_back)
      r = rand()
      if r < mh_rate # accept the proposal
        stat_success_counter += 1
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

# print statistics of the computation

@printf("\nforward_success_counter = %d\nbackward_success_counter = %d\naverage MH rate = %.3f\n", forward_success_counter, backward_success_counter, stat_success_counter * 1.0 / N)

println("\nNo. of solutions in forward rattle")
for i in 1:(max_no_sol+1)
  if stat_num_of_solution_forward[i] > 0
    @printf("%d solutions: %d\n", i-1, stat_num_of_solution_forward[i])
  end
end

println("\nNo. of solutions in backward check")
for i in 1:(max_no_sol+1)
  if stat_num_of_solution_backward[i] > 0
    @printf("%d solutions: %d\n", i-1, stat_num_of_solution_backward[i])
  end
end

# write the sampled states to the file
writedlm("./data/data.txt", sample_data)
