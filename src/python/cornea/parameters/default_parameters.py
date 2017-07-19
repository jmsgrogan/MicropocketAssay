"""
Parameter collections
"""

from copy import deepcopy
from microvessel_chaste.utility import *

import cornea.parameters.parameter_collection
from cornea.parameters.parameter_collection import Parameter 


_default_collection = cornea.parameters.parameter_collection.SimulationParameterCollection()

# Domain Details
_default_collection.add_parameter(Parameter("DomainType", "Planar_2D"))
_default_collection.add_parameter(Parameter("CorneaRadius", 1.3e-3*metres, 
                                           min_val = 1.0, max_val = 1.0, symbol = "R", 
                                           nice_name = "Cornea Radius"))
_default_collection.add_parameter(Parameter("CorneaThickness", 100.0e-6*metres, 
                                           min_val = 1.0, max_val = 1.0, symbol = "T", 
                                           nice_name = "Cornea Thickness"))
_default_collection.add_parameter(Parameter("PelletHeight", 1.0e-3*metres,
                                           min_val = 0.8, max_val = 1.2, symbol = "h", 
                                           nice_name = "Pellet Height"))
_default_collection.add_parameter(Parameter("PelletRadius", 200.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "r", 
                                           nice_name = "Pellet Radius")) 
_default_collection.add_parameter(Parameter("PelletThickness", 40.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "r", 
                                           nice_name = "Pellet Thickness")) 
_default_collection.add_parameter(Parameter("LimbalOffset", 100.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "\epsilon", 
                                           nice_name = "Limbal Offset")) 

# Domain Numerics and Sampling
_default_collection.add_parameter(Parameter("GridSpacing", 40.0e-6*metres, 
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta x", 
                                           nice_name = "Grid Spacing")) 
_default_collection.add_parameter(Parameter("ElementArea2d", 1e3*(1.e-18*metres_cubed)))
_default_collection.add_parameter(Parameter("ElementArea3d", 1e4*(1.e-18*metres_cubed)))                                    
_default_collection.add_parameter(Parameter("NodeSpacing", 40.0e-6*metres, 
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta n", 
                                           nice_name = "Node Spacing"))   
_default_collection.add_parameter(Parameter("DensityGridSpacing", 40.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta x_{\rho}", 
                                           nice_name = "Density Grid Spacing")) 
_default_collection.add_parameter(Parameter("SampleSpacingX", 60.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta x_{s}", 
                                           nice_name = "Sample Spacing X")) 
_default_collection.add_parameter(Parameter("SampleSpacingY", 60.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta y_{s}", 
                                           nice_name = "Sample Spacing Y"))   
_default_collection.add_parameter(Parameter("SampleSpacingZ", 33.0e-6*metres,
                                           min_val = 1.0, max_val = 1.0, symbol = "\Delta z_{s}", 
                                           nice_name = "Sample Spacing Z"))  

# Domain Switches
_default_collection.add_parameter(Parameter("UsePellet", True))  
_default_collection.add_parameter(Parameter("FinitePelletWidth", False)) 

# Angiogenesis details
_default_collection.add_parameter(Parameter("SproutingProbability", (0.5 /3600.0)*per_second, 
                                           min_val = 0.0, max_val = 10.0, symbol = "P_{max}", 
                                           nice_name = "Max Sprouting Probability")) 
_default_collection.add_parameter(Parameter("AttractionStrength", 0.0, 
                                           min_val = 0.0, max_val = 10.0, symbol = "S_{att}", 
                                           nice_name = "Attraction Strength")) 
_default_collection.add_parameter(Parameter("ChemotacticStrength", 0.5, 
                                           min_val = 0.0, max_val = 20.0, symbol = "\chi", 
                                           nice_name = "Chemotactic Strength")) 
_default_collection.add_parameter(Parameter("PersistenceAngle", 5.0, 
                                           min_val = 0.0, max_val = 4.0, symbol = "\alpha", 
                                           nice_name = "Persistence Angle")) 
_default_collection.add_parameter(Parameter("TipVelocity", 20.0 *(1.e-6/3600.0) * metre_per_second, 
                                           min_val = 0.1, max_val = 4.0, symbol = "v", 
                                           nice_name = "Persistence Angle")) 
_default_collection.add_parameter(Parameter("TipExclusionRadius", 40.0e-6*metres, 
                                           min_val = 0.0, max_val = 2.0, symbol = "r_{excl}", 
                                           nice_name = "Tip Exclusion Radius")) 

# Angiogenesis switches
_default_collection.add_parameter(Parameter("DoAnastamosis", True)) 
_default_collection.add_parameter(Parameter("OnlyPerfusedSprout", False)) 

# PDE Details
_default_collection.add_parameter(Parameter("PelletConcentration", 0.3*mole_per_metre_cubed, 
                                           min_val = 0.0, max_val = 1.0, symbol = "C_{0}", 
                                           nice_name = "Initial Pellet Concentration")) 
_default_collection.add_parameter(Parameter("VegfDiffusivity", 6.94e-11 * metre_squared_per_second, 
                                           min_val = 0.1, max_val = 10.0, symbol = "D", 
                                           nice_name = "VEGF Diffusivity")) 
_default_collection.add_parameter(Parameter("VegfDecayRate", (-0.8/3600.0) * per_second, 
                                           min_val = 0.1, max_val = 10.0, symbol = "\lambda", 
                                           nice_name = "VEGF Decay Rate")) 
_default_collection.add_parameter(Parameter("VegfBindingConstant", 100.0, 
                                           min_val = 0.01, max_val = 10.0, symbol = "\theta", 
                                           nice_name = "Pellet Binding Constant")) 
_default_collection.add_parameter(Parameter("VegfBloodConcentration", 0.0*mole_per_metre_cubed, 
                                           min_val = 0.0, max_val = 1.0, symbol = "c_{b}", 
                                           nice_name = "VEGF Blood Concentration")) 
_default_collection.add_parameter(Parameter("VegfPermeability", (3.e-4/3600.0)*metre_per_second, 
                                           min_val = 0.0, max_val = 1.0, symbol = "\rho_{p}", 
                                           nice_name = "Vessel Permeability")) 
_default_collection.add_parameter(Parameter("UptakeRatePerCell", (4.e-18/3600.0)*mole_per_second, 
                                           min_val = 0.0, max_val = 10.0, symbol = "k_{ec}", 
                                           nice_name = "Uptake Rate Per Cell")) 

# PDE Numerics
_default_collection.add_parameter(Parameter("PdeTimeIncrement", 0.01)) 

# PDE Switches
_default_collection.add_parameter(Parameter("IncludeVesselSink", True))
_default_collection.add_parameter(Parameter("UseFixedGradient", False))
_default_collection.add_parameter(Parameter("UsePdeOnly", False))

# Simulation Details
_default_collection.add_parameter(Parameter("TotalTime", 24.0*3600.0*seconds)) 
_default_collection.add_parameter(Parameter("TimeStepSize", 0.5*3600.0*seconds))
_default_collection.add_parameter(Parameter("RunNumber", 0)) 
_default_collection.add_parameter(Parameter("RandomSeed", 0)) 

def get_default_collection():
    return deepcopy(_default_collection)
    