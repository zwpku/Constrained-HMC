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

# solve algebraic equations for given parameters p, with path tracking
function find_solutions_by_tracking(p)
    # Create an empty array.
    S_p = similar(S_p0, 0)
    for s in S_p0
        result = track(tracker, s; target_parameters=p)
        # check that the tracking was successfull
       if is_success(result) && is_real(result)
         sol=solution(result)
	 # check if the solution is new 
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

# find one solution by Newton's method
function find_solution_by_newton(xtmp, grad_xi_vec)
  lam = zeros(k)
  x_now = xtmp
  b = xi(x_now)
  iter = 0
  while norm(b) > newton_res_tol && iter < newton_max_steps
    mat = grad_xi(x_now) * transpose(grad_xi_vec)
    lam += -1.0 * lsmr(mat, b) / step_size 
    x_now = xtmp + step_size * transpose(grad_xi_vec) * lam
    b = xi(x_now)
    iter += 1
  end
  if norm(b) < newton_res_tol
    return reshape(lam, k, 1)
  else 
    return []
  end
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

function find_solutions(p_current, flag)
  if path_tracking_flag == 0
    S_p = find_solutions_total_degreee(p_current)
  else 
    if flag == -1 # prepare the start system, if this is the first time 
      S_p = find_solutions_total_degreee(p_current)
      # record the solutions
      global S_p0 = S_p
      @printf("Starting systems: no. of real solutions = %d\n", length(S_p0))
      #Construct the PathTracker
      global tracker = pathtracker(F; parameters=p, generic_parameters=p_current)
    else 
      S_p = find_solutions_by_tracking(p_current)
    end
  end
  n = length(S_p)
  lambda_vec = zeros(k,n)
  # extract the real part
  for i in 1:n
    lambda_vec[:,i] = [S_p[i][j].re for j in 1:k]
  end
  return lambda_vec
end

