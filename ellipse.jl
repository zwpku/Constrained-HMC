# *********
# This file specifies the model. 
# It defines functions related to ellipse and contains the parameters involved.
# *********

# the reaction coordniate map \xi: \mathbb{R}^d -> \mathbb{R}^k
d = 2
k = 1

# parameters of the ellipse 
c = 3
c2 = c * c

# potential in the target distribution
function V(x)
  return 0.5 * x[1] * x[1]
end

# gradient of the potential V(x)
function grad_V(x)
  return [x[1], 0]
end

# the ith component \xi_i of the map \xi
function xi_i(x, idx)
  return (x[1] * x[1] / c2 + x[2] * x[2]  - 1) * 0.5
end

# gradient of the ith component \xi_i of the map \xi
function grad_xi_i(x, idx)
  return [x[1] / c2 , x[2]]
end

@polyvar lam[1:k] y z p[1:(2 * d)] 

# equations of the Lagrange multipliers, p contains parameters, lam is the unknown multipliers
F = [(p[1] + p[d + 1] * lam[1])^2  / c2 + (p[2] + p[d+2] * lam[1])^2 - 1.0]

#p0 = [3.0,0.0,1.0,1.0]
# parameters for the starting system
p0 = [2.9,-0.6,0.02,-0.18]
