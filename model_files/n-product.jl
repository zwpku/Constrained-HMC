# ************************************
# This file specifies the model. 
# It defines functions related to 3d-torus.
# ************************************

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 2
k = 1

beta = 15.0

# parameters of the ellipse 
c=1.0

# how many different quantities of interest (QoI) will be recorded
num_qoi = 1
# for each quantity of interest, it contains number of bins, lower and upper ranges of the histgram.
qoi_hist_info = [[200, 0.0, 2*pi]] 

# initial state
x0 = ones(d) * c^(1.0/d) 

# user-defined prob. distribution
pj_vec = [[1.0], [0.4, 0.6], [0.2, 0.4, 0.4], [0.2, 0.3, 0.3, 0.2]]

# potential in the target distribution
function V(x)
 return 0.0
end

# gradient of potential V
function grad_V(x)
  return zeros(d)
end

# quantity of interest
function QoI(x)
  return [x[1]]
end

# the ith component \xi_i of the map \xi.
# In this example, \xi is scalar and therefore i=1.

function xi_i(x, idx)
  return prod(x) - c
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  return [c/tmp for tmp in x]
end

if solve_multiple_solutions_frequency > 0 
  # initialize constraint equation for HomotopyContinuation
  @polyvar lam[1:k] p[1:((1+k)*d)]  

  # equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
  #
  # The general form is \xi_i(x + \nabla\xi \lam) = 0, 
  # where \lam is the unknown multipliers, x is given by p[1:d], and \nabla\xi is given by p[(d+1):((k+1)*d)]
  F = [prod([(p[i] + p[d+i]*lam[1]) for i in 1:d]) - c] 
end
