# ************************************
# This specifies the model. 
# The model is the same as the one defined in ellipse.jl, except that it is embeded in R^3 instead of R^2.  
# The purpose of this example is to test the code using a vectorial-valued function \xi$, i.e., k=2.
# ************************************

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 3
k = 2

# parameters of the ellipse 
c = 3
c2 = c * c

# initial state
x0 = [3.0, 0.0, 0.0]

# potential in the target distribution
function V(x)
  return 2.0 * x[1] * x[1]
end

# gradient of the potential V(x)
function grad_V(x)
  return [4.0 * x[1], 0.0, 0.0]
#  return [0.0, 0.0]
end

# the ith component \xi_i of the map \xi
function xi_i(x, idx)
  if idx == 1
    dtmp = (x[1] * x[1] / c2 + x[2] * x[2] + x[3] * x[3] - 1) * 0.5
  else
    dtmp = x[3]
  end
  return dtmp
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  if idx == 1
    return [x[1] / c2 , x[2], x[3]]
  else 
    return [0.0, 0.0, 1.0]
  end
end

if use_homotopy_solver_frequency > 0
  @polyvar lam[1:k] p[1:((1+k)*d)]  

  # equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
  #
  # The general form is \xi_i(x + \nabla\xi \lam) = 0, 
  # where \lam is the unknown multipliers, x is given by p[1:d], and \nabla\xi is given by p[(d+1):((k+1)*d)]

  F = [(p[1] + p[d + 1] * lam[1] + p[2 * d + 1] * lam[2])^2  / c2 + (p[2] + p[d+2] * lam[1] + p[2*d+2] * lam[2])^2 + (p[3] + p[d+3] * lam[1] + p[2*d+3] * lam[2])^2 - 1.0,
	p[3] + p[d+3] * lam[1] + p[2*d+3] * lam[2]]
end
