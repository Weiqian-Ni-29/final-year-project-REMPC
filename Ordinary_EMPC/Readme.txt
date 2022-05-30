This is the realization of Ordinary EMPC in Matlab
The main file is empc.m
make_RPI_set.m generate the RPI set
make_new_disturbances.m generates the predicted disturbances 
noise_generation.m generates disturbances in real simulation
empc_obj.m is the cost function used in fmincon
empc_nonlcon.m is the non-linear constraint used in fmincon
