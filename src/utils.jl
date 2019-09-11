# \xi: a k-dimensional vector
function xi(x)
  return [xi_i(x,i) for i in 1:k]
end

# \nabla\xi: a k\times d matrix
function grad_xi(x, normalize_flag)
  tmp = zeros(k,d)
  for i in 1:k
    tmp[i,:] = grad_xi_i(x,i)
    # normalize the gradient 
    if normalize_flag == 1
      r = norm(tmp[i,:])
      tmp[i,:] /= r
    end
  end
  return tmp
end

# Hamiltonian energy 
function energy(x,v)
  return V(x) + dot(v,v) * 0.5
end

# Solve algebraic equations for given parameters p, with path tracking.
# Note that random number generators are probably NOT used in this function!
function find_solutions_by_tracking(p_current)
    # Create an empty array.
    S_p = similar(S_p0, 0)
    for s in S_p0
       result = track(tracker, s; target_parameters=p_current, accuracy=1e-9)
        # check that the tracking was successfull
       if is_success(result) && is_real(result; tol=1e-8)
         sol=solution(result)
	 # check if the solution is new 
	 new_sol_flag = 1
	 for i in 1:length(S_p)
	    if euclidean_distance(S_p[i], sol) < homotopy_new_sol_tol
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
    mat = grad_xi(x_now, 0) * transpose(grad_xi_vec)
    lam += -1.0 * lsmr(mat, b, atol=newton_matrix_solver_tol, btol=newton_matrix_solver_tol) / step_size 
    x_now = xtmp + step_size * transpose(grad_xi_vec) * lam
    b = xi(x_now)
    iter += 1
  end
  if norm(b) < newton_res_tol
    return reshape(lam, k, 1)
  else # return empty set, no solution has been found
    return Array{Float64}(undef, 0, 0)
  end
end

function find_multiple_solutions(p_current)
  if solve_multiple_solutions_by_homotopy == 1 
    S_p = find_solutions_by_tracking(p_current)
    # check: is the use of length function correct, when k>1?
    n = length(S_p)
    lambda_vec = zeros(k,n)
    # extract the real part
    for i in 1:n
      lambda_vec[:,i] = [S_p[i][j].re for j in 1:k]
    end
  else #instead of using HomotopyContinuation, we use PolynomialRoots package (for k=1)
    # compute the coefficients of polynomial 
    poly = subs(F, p => p_current)
    # there is only one equation, i.e., k=1
    coeff_vec = reverse(poly[1].a)
    # solve the roots
    roots_vec = roots(coeff_vec, epsilon=polyroot_solver_eps)
    sol_vec = []
    for sol in roots_vec
      if abs(sol.im) < 1e-8 
        push!(sol_vec, sol.re)
      end
    end
    lambda_vec = reshape(sol_vec, 1, length(sol_vec))
  end
  return lambda_vec
end

