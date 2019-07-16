import YAML

cfg_file = "./cfg.yml"
println("Reading parameters from file: ", cfg_file)
cfg_data = YAML.load(open(cfg_file))

job_id = cfg_data["job_id"]

# total number of samples 
N = cfg_data["N"]

# the stepsize \tau used in the proposal scheme
step_size = cfg_data["step_size"] 

# tolerance error
check_tol = cfg_data["check_tol"] 

backward_check_tol = cfg_data["backward_check_tol"]

# check the forward scheme and print information, if the flag equals 1.
check_rattle_flag = cfg_data["check_rattle_flag"]

# whether the manifold is algebraic, 
# or, equivalently, whether the function \xi is polynomial 
is_algebraic_manifold_flag = cfg_data["is_algebraic_manifold_flag"]

# how often to solve equation by homotopy method
# Newton's method with be used otherwise
use_homotopy_solver_frequency = cfg_data["use_homotopy_solver_frequency"]

# always use Newton's method, if \xi is not polynomial
if is_algebraic_manifold_flag == 0 && use_homotopy_solver_frequency > 0
  println("Warning: Newton's method will be always used, because manifold is not algebraic!")
  println("\t reset use_homotopy_solver_frequency to zero!")
  use_homotopy_solver_frequency = 0
end

# if homotopy method will be used
if use_homotopy_solver_frequency > 0
  new_sol_tol = cfg_data["new_sol_tol"] 
# the code will be slower, without PathTracking
  path_tracking_flag = cfg_data["path_tracking_flag"]
# upper bound of solution number, here we assume there are at most 4 solutions
end

max_no_sol = cfg_data["max_no_sol"]

# if Newton's solver is used
if use_homotopy_solver_frequency != 1
  # read parameters for Newton's solver
  newton_res_tol = cfg_data["newton_res_tol"]
  newton_max_steps = cfg_data["newton_max_steps"]
end

# if this flag is one, indices are chosen according to the 
# pre-defined probability distributions, based on their distances.
user_defined_pj_flag = cfg_data["user_defined_pj_flag"]

println("Reading parameters from file...done.")
