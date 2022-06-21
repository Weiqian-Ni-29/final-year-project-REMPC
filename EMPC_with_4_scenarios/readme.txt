This is the realization of EMPC with 4 cases in Matlab
The main file is empc4.m
make_RPI_set.m does not play any role in this scenario-based algorithm
make_new_disturbances.m generates the predicted disturbances 
noise_generation.m generates disturbances in real simulation
empc_obj4.m is the cost function used in fmincon
empc_nonlcon4.m is the non-linear constraint used in fmincon
