"""
Classes describing parameters
"""

import os
import random
from copy import deepcopy
import cPickle as pickle
import chaste.core
import microvessel_chaste.utility

class Study():
    
    def __init__(self, work_dir, parameter_collection):
        
        self.work_dir = work_dir
        self.parameter_collection = parameter_collection
        self.range = []
        self.random_realisations = 3
        self.simulation_domains = ["Planar_2D", "Circle_2D", "Planar_3D", "Circle_3D", "Hemisphere_3D"]
        
    def get_task_list(self):
        
        task_list = []
        for eachParameterSet in self.range:
            for idx in range(eachParameterSet[1]):
                if eachParameterSet[1] == 1:
                    param_value = self.parameter_collection.get_parameter(eachParameterSet[0]).value
                else:
                    param = self.parameter_collection.get_parameter(eachParameterSet[0])
                    param_range = (param.max - param.min)*param.value
                    param_value = param.min*param.value + float(idx)/float(eachParameterSet[1])*param_range
                for eachDomain in self.simulation_domains:
                    for jdx in range(self.random_realisations):
                        local_collection = deepcopy(self.parameter_collection)
                        local_collection.get_parameter("domain type").value = eachDomain
                        local_collection.get_parameter("random seed").value = random.randint(0, 1e6)
                        local_collection.get_parameter(eachParameterSet[0]).value = param_value
                        simulation_path = self.work_dir + "/ParamName_" + eachParameterSet[0].replace(" ", "") + "/"
                        simulation_path += "ParamValue_" + str(idx) + "/DomainType_" + eachDomain.replace(" ", "") + "/Run_"+str(jdx)+"/"
                        task_list.append([simulation_path, local_collection])
        return task_list
        
class SimulationParameterCollection:
    
    def __init__(self, random_seed = 1234):
        
        self.collection = {}
        self.random_seed = random_seed
        
    def add_parameter(self, parameter):
        
        self.collection[parameter.name] = parameter
        
    def save(self, file_name):
        with open(file_name, 'wb') as fp:
            pickle.dump([self.collection, self.random_seed], fp)  
            
        # Human friendly version
        output_file = open(os.path.splitext(file_name)[0]+".csv", "w")
        for eachKey in self.collection.keys():
            param = self.collection[eachKey]
            output_file.write(param.name + " , " + str(param.value) + "\n")
        output_file.close()
            
    def load(self, file_name):
        with open(file_name, 'rb') as fp:
            self.collection, self.random_seed = pickle.load(fp) 
            
    def get_parameter(self, name):
        return self.collection[name] 
     
class Parameter:
    
    def __init__(self, name, value, min_val = 1.0, max_val = 1.0, 
                 symbol = None, nice_name = None, lit_source = None):
        
        self.name = name
        self.value = value
        self.min = min_val
        self.max = max_val
        self.value_as_string = ""
        self.store_value_as_string()
        self.unit_dict = {"m" : "metre",
                          "s" : "second",
                          "Hz" : "per_second",
                          "kat" : "mole_per_second",
                          "ms^-1" : "metre_per_second",
                          "m^3" : "metre_cubed",
                          "m^-3mol": "mole_per_metre_cubed",
                          "m^2s^-1" : "metre_squared_per_second"}
        self.symbol = symbol
        self.nice_name = nice_name
        self.lit_source = lit_source
        
    def __getstate__(self):
        self.store_value_as_string()
        d = dict(self.__dict__)
        del d['value']
        return d
    
    def __setstate__(self, d):
        self.__dict__.update(d)
        self.update_value_from_string()
   
    def store_value_as_string(self):
        
        self.value_as_string = str(self.value)
    
    def update_value_from_string(self):
        split_string = self.value_as_string.split()
        if split_string[0] == "False":
            self.value = False
        elif split_string[0] == "True":
            self.value = True 
        elif split_string[0] in ["Planar_2D", "Circle_2D", "Planar_3D", "Circle_3D", "Hemisphere"]:
            self.value = self.value_as_string 
        else:     
            self.value = float(split_string[0])
        if len(split_string)>1 and not split_string[0] in ["Planar_2D", "Circle_2D",  "Planar_3D", "Circle_3D", "Hemisphere"]:
            unit_name = ''.join(split_string[1:])
            self.value*=getattr(microvessel_chaste.utility, self.unit_dict[unit_name])() 
            
if __name__=="__main__":
    
    import cornea.parameters.default_parameters
    
    work_dir = "Python/Cornea/TestParameters/"
    file_handler = chaste.core.OutputFileHandler(work_dir, True)

    collection = cornea.parameters.default_parameters.get_default_collection()    
    collection.save(file_handler.GetOutputDirectoryFullPath()+"/parmaeters.p")
    collection.load(file_handler.GetOutputDirectoryFullPath()+"/parmaeters.p")
    
    for eachParameter in collection.collection.values():
        print eachParameter.value
    