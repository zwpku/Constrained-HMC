import YAML

cfg_file = "./cfg.yml"
println("Reading parameters from file: ", cfg_file)
cfg_data = YAML.load(open(cfg_file))

job_id = cfg_data["job_id"]

# provide filename of the model 
model_file_name = cfg_data["model_file_name"]

# total number of samples 
N = cfg_data["N"]

# output data frequency
output_sample_data_frequency = cfg_data["output_sample_data_frequency"]

# the stepsize \tau used in the proposal scheme
step_size = cfg_data["step_size"] 

# tolerance error
check_tol = cfg_data["check_tol"] 

backward_check_tol = cfg_data["backward_check_tol"]

# check the forward scheme and print information, if the flag equals 1.
check_rattle_flag = cfg_data["check_rattle_flag"]

# the degree of the the constraint equation when it is polynomial.
# it equals zero when it is not polynomial!
degree_polynomial_constraint = cfg_data["degree_polynomial_constraint"]

# how often to solve multiple solutions of the constraint equation 
# Newton's method with be used otherwise
solve_multiple_solutions_frequency = cfg_data["solve_multiple_solutions_frequency"]

# whether use HomotopyContinuation package or not 
solve_multiple_solutions_by_homotopy = cfg_data["solve_multiple_solutions_by_homotopy"]

# always use Newton's method, if the constraint equation is not polynomial
if degree_polynomial_constraint == 0 && solve_multiple_solutions_frequency > 0
  println("Warning: Newton's method will be always used, because the constraint is not polynomial")
  println("\t reset solve_multiple_solutions_frequency to zero!")
  solve_multiple_solutions_frequency = 0
end

@printf("Include model from file: ./%s\n\n", model_file_name)

# assume the model file is under the current directory

include(pwd() * "/" * model_file_name)

if solve_multiple_solutions_frequency > 0
  if k > 1 && solve_multiple_solutions_by_homotopy == 0
    @printf("Warning: constraint equations will be solved by HomotopyContinuation package because k=%d is large than 1.", k)
    println("\t reset solve_multiple_solutions_by_homotopy to 1")
    solve_multiple_solutions_by_homotopy = 1
  end
  # if homotopy method will be used
  if solve_multiple_solutions_by_homotopy == 1
    homotopy_new_sol_tol = cfg_data["homotopy_new_sol_tol"] 
  # the code will be slower, without PathTracking
    path_tracking_in_homotopy_flag = cfg_data["path_tracking_in_homotopy_flag"]
  else # solve scalar constraint equation by PolynomialRoots package.
    polyroot_solver_eps = cfg_data["polyroot_solver_eps"]
  end
end

# upper bound of solution number
max_no_sol = cfg_data["max_no_sol"]
if max_no_sol < degree_polynomial_constraint 
  println("Warning: max_no_sol should be no less than the degree of the polynomial!")
  @printf("\t reset max_no_sol to %d.", degree_polynomial_constraint)
  max_no_sol = degree_polynomial_constraint 
end

# when Newton's solver will be used
if solve_multiple_solutions_frequency != 1
  # read parameters for Newton's solver
  newton_matrix_solver_tol = cfg_data["newton_matrix_solver_tol"]
  newton_res_tol = cfg_data["newton_res_tol"]
  newton_max_steps = cfg_data["newton_max_steps"]
end

# if this flag is one, indices are chosen according to the 
# pre-defined probability distributions, based on their distances.
user_defined_pj_flag = cfg_data["user_defined_pj_flag"]

println("Reading parameters from file...done.")
