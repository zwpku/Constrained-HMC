# ************************************
# This file specifies the model. 
# It defines functions related to 3d-torus.
# ************************************

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 3
k = 1

beta = 20.0

# parameters of the ellipse 
R = 1.0
r = 0.5 

# how many different quantities of interest (QoI) will be recorded
num_qoi = 2
# for each quantity of interest, it contains number of bins, lower and upper ranges of the histgram.
qoi_hist_info = [[200, 0.0, 2*pi], [200, 0.0, 2*pi]] 

# initial state
x0 = [(R+r)/sqrt(2), (R+r)/sqrt(2), 0]

# user-defined prob. distribution
pj_vec = [[1.0], [0.4, 0.6], [0.2, 0.4, 0.4], [0.2, 0.3, 0.3, 0.2]]

# potential in the target distribution
function V(x)
  return (x[1] - x[2])^2 + 5.0 * ((x[1]^2 + x[2]^2)/(R+r)^2 - 1.0)^2 
# return 0.0
end

# gradient of potential V
function grad_V(x)
  return [2.0 * (x[1] - x[2]) + 20.0 * x[1] / (R+r)^2 * ((x[1]^2 + x[2]^2)/(R+r)^2 - 1.0), 2.0 * (x[2] - x[1]) + 20.0 * x[2] / (R+r)^2 * ((x[1]^2 + x[2]^2)/(R+r)^2 - 1.0), 0.0]
end

# quantity of interest
function QoI(x)
  phi = atan(x[3], sqrt(x[1]^2 + x[2]^2) - R) 
  if phi < 0
    phi += 2 * pi
  end
  theta = atan(x[2], x[1])
  if theta < 0
    theta += 2 * pi
  end
  return [phi, theta]
end

# the ith component \xi_i of the map \xi.
# In this example, \xi is scalar and therefore i=1.
# There are different choices of \xi whose zero levelset is 3d-torus. Here, we choose 
# \xi such that it is a fourth order polynomial.
function xi_i(x, idx)
  return (R^2 - r^2 + x[1]^2 + x[2]^2 + x[3]^2)^2 - 4 * R^2 * (x[1]^2 + x[2]^2)
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  tmp = x[1]^2 + x[2]^2 + x[3]^2
  return [4.0 * x[1] * (tmp - r^2 - R^2), 4.0 * x[2] * (tmp - r^2 - R^2), 4.0 * x[3] * (R^2 - r^2 + tmp)]
end

if solve_multiple_solutions_frequency > 0 
  # initialize constraint equation for HomotopyContinuation
  @polyvar lam[1:k] p[1:((1+k)*d)]  

  # equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
  #
  # The general form is \xi_i(x + \nabla\xi \lam) = 0, 
  # where \lam is the unknown multipliers, x is given by p[1:d], and \nabla\xi is given by p[(d+1):((k+1)*d)]
  F = [(R^2 - r^2 + (p[1] + p[d+1] * lam[1])^2 + (p[2] + p[d+2] * lam[1])^2 + (p[3] + p[d+3] * lam[1])^2)^2 - 4 * R^2 * ((p[1] + p[d+1] * lam[1])^2 + (p[2] + p[d+2] * lam[1])^2)]						end
