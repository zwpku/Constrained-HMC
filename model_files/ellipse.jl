# ************************************
# This file specifies the model. 
# It defines functions related to ellipse and contains the parameters in the system.
# ************************************

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 2
k = 1

# parameters of the ellipse 
c = 3
c2 = c * c

# initial state
x0 = [3.0, 0.0]

# potential in the target distribution
function V(x)
#  return 2.0 * x[1] * x[1]
  return 0.0 
end

# gradient of the potential V(x)
function grad_V(x)
#  return [4.0 * x[1], 0.0]
  return [0.0, 0.0]
end

# the ith component \xi_i of the map \xi
function xi_i(x, idx)
  return (x[1] * x[1] / c2 + x[2] * x[2]  - 1.0) * 0.5
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  return [x[1] / c2 , x[2]]
end

if use_homotopy_solver_frequency > 0
  @polyvar lam[1:k] p[1:((1+k)*d)]  

  # equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
  #
  # The general form is \xi_i(x + \nabla\xi \lam) = 0, 
  # where \lam is the unknown multipliers, x is given by p[1:d], and \nabla\xi is given by p[(d+1):((k+1)*d)]
  F = [(p[1] + p[d + 1] * lam[1])^2  / c2 + (p[2] + p[d+2] * lam[1])^2 - 1.0]
end
