# ************************************
# This file specifies the model. 
# ************************************

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 10
k = 2

beta = 1.0

# parameters of the ellipse 
R = 3.0
c= 2.0

# how many different quantities of interest (QoI) will be recorded
num_qoi = 5
# for each quantity of interest, it contains number of bins, lower and upper ranges of the histgram.
qoi_hist_info = [[4, -0.5, 3.5], [300, -3.0, 3.0], [300, -3.0, 3.0], [300, -3.0, 3.0], [300, -3.0, 3.0]] 

# initial state
x0 = ones(d)  
x0[1] = 2.0

# user-defined prob. distribution
pj_vec = [[1.0], [0.4, 0.6], [0.2, 0.4, 0.4], [0.2, 0.3, 0.3, 0.2]]

# potential in the target distribution
function V(x)
 return 0.5 * (x[1]-0.6)^2
end

# gradient of potential V
function grad_V(x)
  tmp = zeros(d)
  tmp[1] = x[1] - 0.6
  return tmp
end

# quantity of interest
function QoI(x)
  idx = 1 * (x[1] > 0) + 2 * (x[2] > 0) + 3 * (x[3] > 0)
  if idx == 6
    idx = 0
  end
  return [idx, x[1], x[2], x[3], x[4]]
end

# the ith component \xi_i of the map \xi.
# In this example, k=2, therefore, i=1 or 2.

function xi_i(x, idx)
  if idx == 1
    return prod(x[1:3]) - c
  else
    return 0.5 * (sum(x .^ 2) - R^2)
  end
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  if idx == 1
    tmp_vec = zeros(d)
    tot_prod = prod(x[1:3])
    tmp_vec[1:3] = [tot_prod/x[i] for i in 1:3]
    return tmp_vec
  else 
    return x
  end
end

if solve_multiple_solutions_frequency > 0 
  # initialize constraint equation for HomotopyContinuation
  @polyvar lam[1:k] p[1:((1+k)*d)]  

  # equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
  #
  # The general form is \xi_i(x + \nabla\xi \lam) = 0, 
  # where \lam is the unknown multipliers, x is given by p[1:d], and \nabla\xi is given by p[(d+1):((k+1)*d)]

  F = [prod([(p[i] + p[d+i]*lam[1] + p[2*d+i]*lam[2]) for i in 1:3]) - c, 0.5 * (sum([(p[i] + p[d+i]*lam[1] + p[2*d+i]*lam[2])^2 for i in 1:d]) - R^2) ] 
end
