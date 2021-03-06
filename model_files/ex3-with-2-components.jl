# ************************************
# This file specifies the model. 
# ************************************

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 5
k = 1

beta = 1.0

# parameters of the ellipse 
c=2.0

# how many different quantities of interest (QoI) will be recorded
num_qoi = 1
# for each quantity of interest, it contains number of bins, lower and upper ranges of the histgram.
qoi_hist_info = [[200, 0.0, 2*pi]] 

# initial state
x0 = zeros(d)  
x0[1] = c
x0[d] = 1

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
  return 0.5 * ((x[1]^2 - c^2)^2 + sum(x[2:d] .^ 2) - 1.0)
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  tmp = zeros(d)
  tmp[1] = 2 * x[1] * (x[1]^2 - c^2)
  tmp[2:d] = x[2:d]
  return tmp
end

if solve_multiple_solutions_frequency > 0 
  # initialize constraint equation for HomotopyContinuation
  @polyvar lam[1:k] p[1:((1+k)*d)]  

  # equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
  #
  # The general form is \xi_i(x + \nabla\xi \lam) = 0, 
  # where \lam is the unknown multipliers, x is given by p[1:d], and \nabla\xi is given by p[(d+1):((k+1)*d)]

  F = [0.5 * (((p[1] + p[d+1]*lam[1])^2 - c^2)^2 + sum([(p[i] + p[d+i]*lam[1])^2 for i in 2:d]) - 1.0)] 
end
