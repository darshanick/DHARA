from compressible import Compressible
from fns import time_advance_euler, time_advance_RK2
import time 
from boundary import boundary

ti = time.time()                    # Initial time

compress = Compressible()           # Making of Object
compress.set_arrays()               # Setting of variables
compress.init_hydro()               # Initiating inital profiles at time tinit
boundary(compress)                  # Boundary conditions on initial primitive variables 
compress.density()                  # Finding initial rho using initial p & T after proper boundary
compress.specific_energy()          # Finding initial e_t using initial ux, uz & T after proper boundary  

# print ('ti = ', compress.rho)

# time_advance_euler(compress)        # Time advance steps using Euler method
time_advance_RK2(compress)          # Time advance steps using RK2 method

# print ('tf = ', compress.rho)

tf = time.time()                    # Final time

print('Total time of simulation = ' ,tf-ti)          # Time taken to run the code


