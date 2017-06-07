"""
Parameter collections
"""

from copy import deepcopy
from microvessel_chaste.utility import *

import cornea.parameters.parameter_collection
from cornea.parameters.parameter_collection import Parameter 


_default_collection = cornea.parameters.parameter_collection.SimulationParameterCollection()

# Domain Details
_default_collection.add_parameter(Parameter("domain type", "Planar 2D"))
_default_collection.add_parameter(Parameter("cornea radius", 1.3e-3*metre(), 
                                           min_val = 1.0, max_val = 1.0, symbol = "R", 
                                           nice_name = "Cornea Radius"))
_default_collection.add_parameter(Parameter("cornea thickness", 100.0e-6*metre(), 
                                           min_val = 1.0, max_val = 1.0, symbol = "T", 
                                           nice_name = "Cornea Thickness"))
_default_collection.add_parameter(Parameter("pellet height", 1.0e-3*metre(),
                                           min_val = 0.8, max_val = 1.2, symbol = "h", 
                                           nice_name = "Pellet Height"))
_default_collection.add_parameter(Parameter("pellet radius", 300.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "r", 
                                           nice_name = "Pellet Radius")) 
_default_collection.add_parameter(Parameter("pellet thickness", 40.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "r", 
                                           nice_name = "Pellet Thickness")) 
_default_collection.add_parameter(Parameter("limbal offset", 200.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "\epsilon", 
                                           nice_name = "Limbal Offset")) 

# Domain Numerics and Sampling
_default_collection.add_parameter(Parameter("grid spacing", 40.0e-6*metre(), 
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta x", 
                                           nice_name = "Grid Spacing")) 
_default_collection.add_parameter(Parameter("element area 2d", 1e3*(1.e-18*metre_cubed())))
_default_collection.add_parameter(Parameter("element area 3d", 1e4*(1.e-18*metre_cubed())))                                    
_default_collection.add_parameter(Parameter("node spacing", 40.0e-6*metre(), 
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta n", 
                                           nice_name = "Node Spacing"))   
_default_collection.add_parameter(Parameter("density grid spacing", 40.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta x_{\rho}", 
                                           nice_name = "Density Grid Spacing")) 
_default_collection.add_parameter(Parameter("sample spacing x", 20.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta x_{s}", 
                                           nice_name = "Sample Spacing X")) 
_default_collection.add_parameter(Parameter("sample spacing y", 20.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta y_{s}", 
                                           nice_name = "Sample Spacing Y"))   
_default_collection.add_parameter(Parameter("sample spacing z", 20.0e-6*metre(),
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta z_{s}", 
                                           nice_name = "Sample Spacing Z"))  

# Domain Switches
_default_collection.add_parameter(Parameter("use pellet", False))  
_default_collection.add_parameter(Parameter("use finite pellet width", False)) 

# Angiogenesis details
_default_collection.add_parameter(Parameter("sprouting probability", 0.5 /(3600.0*second()), 
                                           min_val = 0.0, max_val = 10.0, symbol = "P_{max}", 
                                           nice_name = "Max Sprouting Probability")) 
_default_collection.add_parameter(Parameter("attraction strength", 0.0, 
                                           min_val = 0.0, max_val = 10.0, symbol = "S_{att}", 
                                           nice_name = "Attraction Strength")) 
_default_collection.add_parameter(Parameter("chemotactic strength", 0.5, 
                                           min_val = 0.0, max_val = 20.0, symbol = "\chi", 
                                           nice_name = "Chemotactic Strength")) 
_default_collection.add_parameter(Parameter("persistence angle", 5.0, 
                                           min_val = 0.0, max_val = 4.0, symbol = "\alpha", 
                                           nice_name = "Persistence Angle")) 
_default_collection.add_parameter(Parameter("tip exclusion radius", 40.0e-6*metre(), 
                                           min_val = 0.0, max_val = 2.0, symbol = "r_{excl}", 
                                           nice_name = "Tip Exclusion Radius")) 

# Angiogenesis switches
_default_collection.add_parameter(Parameter("do anastomosis", True)) 

# PDE Details
_default_collection.add_parameter(Parameter("pellet concentration", 0.3*mole_per_metre_cubed(), 
                                           min_val = 0.0, max_val = 1.0, symbol = "C_{0}", 
                                           nice_name = "Initial Pellet Concentration")) 
_default_collection.add_parameter(Parameter("vegf diffusivity", 6.94e-11 * metre_squared_per_second(), 
                                           min_val = 0.1, max_val = 10.0, symbol = "D", 
                                           nice_name = "VEGF Diffusivity")) 
_default_collection.add_parameter(Parameter("vegf decay rate", (-0.8/3600.0) * per_second(), 
                                           min_val = 0.1, max_val = 10.0, symbol = "\lambda", 
                                           nice_name = "VEGF Decay Rate")) 
_default_collection.add_parameter(Parameter("pellet binding constant", 100.0, 
                                           min_val = 0.01, max_val = 10.0, symbol = "\theta", 
                                           nice_name = "Pellet Binding Constant")) 
_default_collection.add_parameter(Parameter("vegf blood concentration", 0.0*mole_per_metre_cubed(), 
                                           min_val = 0.0, max_val = 1.0, symbol = "c_{b}", 
                                           nice_name = "VEGF Blood Concentration")) 
_default_collection.add_parameter(Parameter("vessel permeability", (3.e-4/3600.0)*metre_per_second(), 
                                           min_val = 0.0, max_val = 1.0, symbol = "\rho_{p}", 
                                           nice_name = "Vessel Permeability")) 
_default_collection.add_parameter(Parameter("uptake rate per cell", (5.e8/3600.0)*mole_per_second(), 
                                           min_val = 0.0, max_val = 10.0, symbol = "k_{ec}", 
                                           nice_name = "Uptake Rate Per Cell")) 

# PDE Numerics
_default_collection.add_parameter(Parameter("pde time increment", 0.01)) 

# PDE Switches
_default_collection.add_parameter(Parameter("include vessel sink", True))
_default_collection.add_parameter(Parameter("use fixed gradient", False))
_default_collection.add_parameter(Parameter("use pde only", False))

# Simulation Details
_default_collection.add_parameter(Parameter("total time", 24)) 
_default_collection.add_parameter(Parameter("time step", 0.5)) 
_default_collection.add_parameter(Parameter("run number", 0)) 
_default_collection.add_parameter(Parameter("random seed", 0)) 

def get_default_collection():
    return deepcopy(_default_collection)
    