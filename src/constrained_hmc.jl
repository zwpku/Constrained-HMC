using HomotopyContinuation, DynamicPolynomials
using PolynomialRoots
using LinearAlgebra
using IterativeSolvers 
using PyPlot
using DelimitedFiles
using Printf

include("read_params.jl")

include("utils.jl")

# compute several possible proposal states, at the current state (x,v)
function forward_rattle(x, v, use_newton_flag)
  grad_pot_vec = grad_V(x)
  # this should be a (k x d) matrix
  grad_xi_vec = grad_xi(x)
  coeff = - step_size * step_size * 0.5
  x_tmp = x + step_size * v + coeff * grad_pot_vec 
  # find Lagrange multipliers for x
  if use_newton_flag == 1 # by Newton's method
    lam_x = find_solution_by_newton(x_tmp, grad_xi_vec)
  else  # find multiple solutions 
    # be careful how the parameters are ordered in p
    p = vcat(x_tmp, step_size * reshape(transpose(grad_xi_vec), length(grad_xi_vec), 1)[:,1])
    lam_x = find_multiple_solutions(p, use_newton_flag)
  end
  n = size(lam_x, 2)
  # if we find at least one solutions
  if n > 0
    if n == 1 # if there is only one solution
      j = 1
      pj = 1.0
    else # when multiple solutions exist
      if user_defined_pj_flag == 1 
	# sort by distances, and choose indices according to the probability distribution in pj_vec.
	dist = zeros(n)
	# compute the distance wrt x for each proposal
	for i in 1:n
	  dist[i] = norm(x_tmp + step_size * transpose(grad_xi_vec) * lam_x[:,i] - x)
	end
	# sort accrording to the distance
	perm = sortperm(dist)
	# generate a random number uniformly in [0,1]
	r = rand()
	# decide which index to use
	jtmp = searchsortedfirst(pj_acc_vec[n], r)
	j = perm[jtmp]
	# pj is the corresponding probability
	pj = pj_vec[n][jtmp]
      else 
	 # choose one index uniformly
	j = rand(1:n)
	pj = 1.0 / n
      end
    end
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
    # Note: matrix inversion may be numerically problematic, when k is large! Consider replacing it by iterative solver, such as lsmr.
    lam_v = - inv(mat_v_tmp) * grad_xi_vec_1 * v_tmp 
    # compute the updated velocity v^1
    v_1 = v_tmp + transpose(grad_xi_vec_1) * lam_v
    return n, pj, x_1, v_1
  end # no solutions are found, if we reach here
  return 0, 0, x, v
end

#backward check, similar to the forward_rattle function.
function backward_check(x1, v1, x, v, use_newton_flag)
  grad_pot_vec = grad_V(x1)
  # this should be a (k x d) matrix
  grad_xi_vec = grad_xi(x1)
  coeff = - step_size * step_size * 0.5
  x_tmp = x1 + step_size * v1 + coeff * grad_pot_vec 
  if use_newton_flag == 1 # by Newton's method
    lam_x = find_solution_by_newton(x_tmp, grad_xi_vec)
  else  # find multiple solutions
    # be careful how the parameters are ordered in p
    p = vcat(x_tmp, step_size * reshape(transpose(grad_xi_vec), length(grad_xi_vec), 1)[:,1]) 
    lam_x = find_multiple_solutions(p, use_newton_flag)
  end
  n_back = size(lam_x, 2)
  backward_found_flag = 0
  pj_back = 0.0
  # sort the solutions by distance
  if n_back > 1 && user_defined_pj_flag == 1 
    # sort by distances
    dist = zeros(n_back)
    # compute the distance wrt x1 for each proposal
    for i in 1:n_back
      dist[i] = norm(x_tmp + step_size * transpose(grad_xi_vec) * lam_x[:,i] - x1)
    end
    # sort accrording to the distance
    perm = sortperm(dist)
  else 
    perm = 1:n_back
  end
  # go through all solutions, and check if one solution is (x,v)
  for j in 1:n_back
    # compute the new state x^2
    x_2 = x_tmp + step_size * transpose(grad_xi_vec) * lam_x[:,perm[j]]
    # check whether the states are the same
    # Note that, we don't check velocities, because the velocities will be the same, if the states are the same.   
    if norm(x_2 - x) < backward_check_tol 
      backward_found_flag = 1
      # record the probability as well. It will be used to compute the M-H rate
      if n_back > 1
	pj_back = pj_vec[n_back][j]
      else 
        pj_back = 1.0
      end
      break 
    end
  end
  return backward_found_flag, n_back, pj_back
end

function rand_draw_velocity(x)
  grad_xi_vec = grad_xi(x)
  # generate the orthnormal basis of the tangent space
  U_x = nullspace(grad_xi_vec)
  # generate normal Gaussian vector, as coefficients under the basis  
  coeff = randn(d-k)
  return U_x * coeff
end

# array to store the samples 
sample_data_vec = []
forward_success_counter = 0
backward_success_counter = 0
newton_counter = 0 
stat_success_counter = 0
stat_average_distance = 0

stat_num_of_solution_forward = zeros(max_no_sol+1)
stat_num_of_solution_backward = zeros(max_no_sol+1)

# when mutilple solutions will be solved 
if solve_multiple_solutions_frequency > 0
  # initialize the vector pj_vec 
  if user_defined_pj_flag == 1
    pj_vec = [[1.0], [0.0, 1.0], [0.6, 0.3, 0.1], [0.6, 0.2, 0.1, 0.1]]
  else #uniform distribution
    pj_vec = [[1.0 / i for j in 1:i] for i in 1:max_no_sol]
  end
  pj_acc_vec = similar(pj_vec)
  for i in 1:max_no_sol
    pj_acc_vec[i] = cumsum(pj_vec[i])
  end

  if solve_multiple_solutions_by_homotopy == 1  
    # prepare the start system
    global num_sol_start_system = 0
    # find a system as many solutions as possible
    for i in 1:10
      v0 = rand_draw_velocity(x0)
      # -1 indicates that we are solving the start system
      use_newton_flag = -1
      n, pj, x1, v1 = forward_rattle(x0, v0, use_newton_flag)
      if n > 0
	global x0 = x1
      end
    end
    if num_sol_start_system > 0
      @printf("Starting systems: no. of solutions = %d\n", length(S_p0))
      println("Starting systems: x0=", x0)
    else 
      println("Error: couldn't find a start system with nonzero solutions!")
      exit(1)
    end
  end
end

if num_qoi <= 0
  @printf("\nWarning: no quantity of interest will be computed (num_qoi = %d)!\n", num_qoi)
else 
  qoi_data_vec = []
  # record the bin width 
  for i in 1:num_qoi
    push!(qoi_hist_info[i], (qoi_hist_info[i][3] - qoi_hist_info[i][2]) / qoi_hist_info[i][1])
  end
  qoi_counter = [zeros(trunc(Int32, qoi_hist_info[i][1])) for i in 1:num_qoi]
end

@printf("\nSampling started... \n")
flush(stdout)

# count the runtime 
start_time = time()
@time begin
# the main loop 
for i in 1:N
  global x0, v0, forward_success_counter, backward_success_counter, stat_success_counter, stat_num_of_solution_forward, stat_num_of_solution_backward, stat_average_distance

  # output some information 
  if N > 100 && i % (div(N, 10)) == 0
    current_time = time()
    @printf(" Step %d, %.1f%% finished.", i, i * 100 / N)
    @printf("\n\tElapsed time = %.2fSec., (Estimated) Total time = %.2fSec., Remaining time = %.2fSec.", (current_time - start_time), (current_time - start_time) / i * N, (current_time - start_time) / i * (N-i))
    @printf("\n\tForward_success_counter = %d (%.1f%%)\n\tBackward_success_counter = %d (%.1f%%)\n", forward_success_counter, forward_success_counter * 100 / i, backward_success_counter, backward_success_counter * 100 / forward_success_counter)
    flush(stdout)
  end 
  # first of all, randomly update the velocity 
  v0 = rand_draw_velocity(x0)

  if output_sample_data_frequency > 0 && i % output_sample_data_frequency == 0 # record the current state for output
    push!(sample_data_vec, vcat(x0, v0))
  end

  if num_qoi > 0
    qoi_val = QoI(x0)

    if output_qoi_data_frequency > 0 && i % output_qoi_data_frequency == 0 # record the current QoI
      push!(qoi_data_vec, qoi_val)
    end

    for j in 1:num_qoi
      idx = trunc(Int32, (qoi_val[j] - qoi_hist_info[j][2]) / qoi_hist_info[j][4] ) + 1
      if idx < 1 
        idx = 1
      end 
      if idx > qoi_hist_info[j][1] 
        idx = qoi_hist_info[j][1]
      end 
      qoi_counter[j][idx] += 1
    end
  end

  if solve_multiple_solutions_frequency > 0 && i % solve_multiple_solutions_frequency == 0
    use_newton_flag = 0
  else 
    use_newton_flag = 1
    global newton_counter += 1
  end
  # compute proposal states
  n, pj, x1, v1 = forward_rattle(x0, v0, use_newton_flag)
  if n <= max_no_sol 
    stat_num_of_solution_forward[n+1] += 1
  else 
    @printf("Warning: No. of solutions in forward rattle (=%d) is larger than upper bound (=%d)!", n, max_no_sol)
  end
  if n == 0 # no solution
    continue
  end
  # otherwise, if one solution has been found, do backward check
  forward_success_counter += 1 
  # reverse the velocity, and do backward check
  found_flag, n_back, pj_back = backward_check(x1, -v1, x0, -v0, use_newton_flag)
  if n_back <= max_no_sol 
    stat_num_of_solution_backward[n_back+1] += 1
  else 
    @printf("Warning: No. of solutions in backward check (=%d) is larger than upper bound (=%d)!", n_back, max_no_sol)
  end
  if found_flag == 1 # if the backward check is passed 
    backward_success_counter += 1
    h = energy(x0, v0) 
    h_1 = energy(x1, v1) 
#      @printf("n=%d, pj = %.3f n_back=%d, pj_back=%.3f\n", n, pj, n_back, pj_back)
    # compute the MH-rate
    mh_rate = min(1.0, exp(h - h_1) * pj_back / pj)
    r = rand()
    if r < mh_rate # accept the proposal
      stat_success_counter += 1
      stat_average_distance += norm(x1-x0)
      x0 = x1
      v0 = v1
      continue
    end
  end
# it is not necessary to flip the velocity, 
# because the velocity will be redrawn at the begining of the next step
#  v0 = -v0
end
# time end
end

@printf("Sampling finished.\n")
# print statistics of the computation

@printf("\n========== Statistics ==========")

@printf("\nForward_success_counter = %d (%.1f%%)\nBackward_success_counter = %d (%.1f%%)\n", forward_success_counter, forward_success_counter * 100 / N, backward_success_counter, backward_success_counter * 100 / forward_success_counter)

@printf("Average MH rate = %.3f\nTotal successful jump rate = %.3f\nAverage jump distance (including no jump)= %.3f (%.3f)\n", stat_success_counter * 1.0 / backward_success_counter, stat_success_counter / N, stat_average_distance * 1.0 / stat_success_counter, stat_average_distance * 1.0 / N)

@printf("\nNo. of steps using Newton's method: %d\n", newton_counter)

if solve_multiple_solutions_by_homotopy == 1
  @printf("No. of steps using HomotopyContinuation package: %d\n", N - newton_counter)
else 
  @printf("No. of steps using PolynomialRoots package: %d\n", N - newton_counter)
end

println("\nNo. of solutions in forward rattle:")
for i in 1:(max_no_sol+1)
  if stat_num_of_solution_forward[i] > 0
    @printf("\t%d solutions: %d (%.1f%%)\n", i-1, stat_num_of_solution_forward[i], stat_num_of_solution_forward[i] * 100 / N)
  end
end

println("\nNo. of solutions in backward check:")
for i in 1:(max_no_sol+1)
  if stat_num_of_solution_backward[i] > 0
    @printf("\t%d solutions: %d (%.1f%%)\n", i-1, stat_num_of_solution_backward[i], stat_num_of_solution_backward[i] * 100 / forward_success_counter)
  end
end

if output_sample_data_frequency > 0
  # write the recorded states to file
  writedlm(@sprintf("./traj_data_%d.txt", job_id), sample_data_vec)
end

if num_qoi > 0
  # write the quantity of interest (QoI) to file
  if output_qoi_data_frequency > 0
    writedlm(@sprintf("./qoi_data_%d.txt", job_id), qoi_data_vec)
  end
  # write the hist of quantity of interest (QoI) to file
  output_file = open(@sprintf("./qoi_counter_%d.txt", job_id), "w")
  writedlm(output_file, [N num_qoi])
  writedlm(output_file, qoi_hist_info)
  writedlm(output_file, qoi_counter)
  close(output_file)
end