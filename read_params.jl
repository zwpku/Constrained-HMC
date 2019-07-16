import YAML

cfg_file = "./cfg.yml"
println("Reading parameters from file: ", cfg_file)
cfg_data = YAML.load(open(cfg_file))

job_id = cfg_data["job_id"]

# the stepsize \tau used in the proposal scheme
step_size = cfg_data["step_size"] 

# check the forward scheme and print information, if the flag equals 1.
check_rattle_flag = cfg_data["check_rattle_flag"]

# tolerance error
check_tol = cfg_data["check_tol"] 

backward_check_tol = cfg_data["backward_check_tol"]
new_sol_tol = cfg_data["new_sol_tol"] 

# newton solver
newton_res_tol = cfg_data["newton_res_tol"]
newton_max_steps = cfg_data["newton_max_steps"]

# how often to solve equation by homotopy method
# Newton's method with be used otherwise
use_homotopy_solver_frequency = cfg_data["use_homotopy_solver_frequency"]

# the code will be slower, without PathTracking
path_tracking_flag = cfg_data["path_tracking_flag"]

# total number of samples 
N = cfg_data["N"]

# upper bound of solution number, here we assume there are at most 4 solutions
max_no_sol = cfg_data["max_no_sol"]

# if this flag is one, indices are chosen according to the 
# pre-defined probability distributions, based on their distances.
user_defined_pj_flag = cfg_data["user_defined_pj_flag"]

println("Reading parameters from file...done.")
