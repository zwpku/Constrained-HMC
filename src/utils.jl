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

function ode_f(u,p,t)
  return -1.0 * transpose(grad_xi(u, 0)) * xi(u)
end

function find_initial_state_by_ODE(x0)
  tspan = (0.0, 100.0)
  prob = ODEProblem(ode_f, x0, tspan)
  sol = OrdinaryDiffEq.solve(prob, DynamicSS(Tsit5()), save_everystep=false)
  x0 = sol[2]
  return x0
end

# Solve algebraic equations for given parameters p, with path tracking.
# Note that random number generators are probably NOT used in this function!
function find_solutions_by_tracking(p_current)
    # Create an empty array.
    Gp = [subs(f, p => p_current) for f in F]
    R = HomotopyContinuation.solve(F_p, Gp, S_p0)
    S_p = solutions(R, only_real=true, only_nonsingular=true)
#    println("result=", S_p)
#=
    S_p = similar(S_p0, 0)
    for s in S_p0
       result = track(tracker, s; target_parameters=p_current, accuracy=1e-8)
        # check that the tracking was successfull
       if issuccess(result) 
	 push!(S_p, solution(result))
       end
    end
=#
    return S_p
end

# find one solution by Newton's method
function find_solution_by_newton(xtmp, grad_xi_vec; res_tol=newton_res_tol, max_steps=newton_max_steps, lam0=zeros(k))
  lam=lam0
  x_now = xtmp + step_size * transpose(grad_xi_vec) * lam
  b = xi(x_now)
  iter = 0
  while norm(b) > res_tol && iter < max_steps
    mat = grad_xi(x_now, 0) * transpose(grad_xi_vec)
    lam += -1.0 * lsmr(mat, b, atol=newton_matrix_solver_tol, btol=newton_matrix_solver_tol) / step_size
    x_now = xtmp + step_size * transpose(grad_xi_vec) * lam
    b = xi(x_now)
    iter += 1
  end
  if norm(b) < res_tol
    return iter, reshape(lam, k, 1) 
  else # return empty set, no solution has been found
    return iter, Array{Float64}(undef, 0, 0)
  end
end

function find_multiple_solutions(x_tmp, grad_xi_vec)
  p_current = vcat(x_tmp, reshape(transpose(grad_xi_vec), length(grad_xi_vec), 1)[:,1])
  if solve_multiple_solutions_by_homotopy == 1 
    roots_vec = find_solutions_by_tracking(p_current)
  else #instead of using HomotopyContinuation, we use PolynomialRoots package (for k=1)
    # compute the coefficients of polynomial 
    poly = subs(F, p => p_current)
    # there is only one equation, i.e., k=1
    coeff_vec = reverse(poly[1].a)
    # solve the roots
    roots_vec = roots(coeff_vec, epsilon=polyroot_solver_eps)
  end
  n = length(roots_vec)
  prev_sol_vec = []
  for idx in 1:n
    tmp = norm([roots_vec[idx][i].im for i in 1:k])
    # only choose real solutions
    if tmp > 1e-8 
      continue
    end
    tmp_real = [roots_vec[idx][i].re / step_size for i in 1:k]
    # refine the solution using Newton's method
    if refine_by_newton_max_step > 0
      num, tmp_real = find_solution_by_newton(x_tmp, grad_xi_vec, res_tol= check_tol, max_steps=refine_by_newton_max_step, lam0=tmp_real)
      global stat_tot_refine_newton_steps += num 
    end
    if length(tmp_real) == k
      # only if the current solution is a new one
      new_flag = 1
      for sol in prev_sol_vec
	if norm(tmp_real - sol) < new_sol_tol
	   new_flag = 0
	   break
	end
      end
      if new_flag == 1 
	push!(prev_sol_vec, tmp_real)
      end
    end
  end
  n = length(prev_sol_vec)
  sol_vec = zeros(k, n)
  for ii in 1:n
    sol_vec[:, ii] = prev_sol_vec[ii]
  end
  return sol_vec 
end

