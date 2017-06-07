
"""Copyright (c) 2005-2017, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import chaste # Core Chaste functionality
chaste.init() # Initialize MPI and PETSc

import cornea.simulations.base_simulation
import cornea.parameters.default_parameters
        
if __name__ == '__main__':
    
    work_dir = "Python/Cornea/TestSimulationPde/"
    parameter_collection = cornea.parameters.default_parameters.get_default_collection()
    parameter_collection.get_parameter("include vessel sink").value = False
    
    domain_types = ["Planar 2D", "Circle 2D", "Planar 3D", "Circle 3D", "Hemisphere 3D"]
    #domain_types = ["Planar 2D"]
    random_seeds = [1234]

    for eachDomainType in domain_types:
        run_number = 0
        parameter_collection.get_parameter("domain type").value = eachDomainType
        for eachSeed in random_seeds:
            parameter_collection.get_parameter("run number").value = run_number
            parameter_collection.get_parameter("random seed").value = eachSeed
            local_work_dir = work_dir + "/DomainType_" + eachDomainType.replace(" ", "") + "/Run_" + str(run_number)
            simulation = cornea.simulations.base_simulation.BaseSimulation(parameter_collection, work_dir)
            simulation.run()
            run_number += 1
    