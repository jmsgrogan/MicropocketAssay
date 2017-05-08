"""
Classes describing parameters
"""

import cPickle as pickle
import chaste.core
import microvessel_chaste.utility

class SimulationParameterCollection:
    
    def __init__(self, random_seed = 1234):
        
        self.collection = {}
        self.random_seed = random_seed
        
    def add_parameter(self, parameter):
        
        self.collection[parameter.name] = parameter
        
    def save(self, file_name):
        with open(file_name, 'wb') as fp:
            pickle.dump([self.collection, self.random_seed], fp)  
            
    def load(self, file_name):
        with open(file_name, 'rb') as fp:
            self.collection, self.random_seed = pickle.load(fp) 
            
    def get_parameter(self, name):
        return collection[name] 
     
class Parameter:
    
    def __init__(self, name, value, min_val = None, 
                 max_val = None, symbol = None, nice_name = None):
        
        self.name = name
        self.value = value
        self.min = min_val
        self.max = max_val
        self.value_as_string = ""
        self.store_value_as_string()
        self.unit_dict = {"m" : "metre"}
        self.symbol = symbol
        self.nice_name = nice_name
        
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
        else:           
            self.value = float(split_string[0])
        if len(split_string)==2:
            self.value*=getattr(microvessel_chaste.utility, self.unit_dict[split_string[1]])() 
            
if __name__=="__main__":
    
    from microvessel_chaste.utility import *
    
    work_dir = "Python/Cornea/TestParameters/"
    file_handler = chaste.core.OutputFileHandler(work_dir, True)

    collection = SimulationParameterCollection()
    collection.add_parameter(Parameter("pellet height", 1.0e-3*metre()))
    collection.add_parameter(Parameter("cornea radius", 1.3e-3*metre()))
    collection.add_parameter(Parameter("cornea thickness", 80.0e-6*metre()))    
    collection.add_parameter(Parameter("pellet thickness", 40.0e-6*metre()))
    collection.add_parameter(Parameter("grid spacing", 40.0e-6*metre()))  
    collection.add_parameter(Parameter("node spacing", 40.0e-6*metre()))  
    collection.add_parameter(Parameter("limbal offset", 200.0e-6*metre()))  
    collection.add_parameter(Parameter("density grid spacing", 40.0e-6*metre())) 
    collection.add_parameter(Parameter("sample spacing x", 20.0e-6*metre()))  
    collection.add_parameter(Parameter("sample spacing y", 20.0e-6*metre()))  
    collection.add_parameter(Parameter("sample spacing z", 20.0e-6*metre()))  
    collection.add_parameter(Parameter("use pellet", False))  
    collection.add_parameter(Parameter("pellet radius", 300.0e-6*metre()))  
    collection.add_parameter(Parameter("use finite pellet width", False))          
    
    collection.save(file_handler.GetOutputDirectoryFullPath()+"/parmaeters.p")
    collection.load(file_handler.GetOutputDirectoryFullPath()+"/parmaeters.p")
    
    for eachParameter in collection.collection:
        print eachParameter.value
    