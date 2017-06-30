import numpy as np
import matplotlib.pyplot as plt
from symfit import parameters, variables, Fit, exp, Parameter
import vtk
import chaste.core

import cornea.postprocessing.sampled_quantity

def plot_3d(work_dir, coords, vals, name="Default"):

    polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    point_data = vtk.vtkDoubleArray()
    point_data.SetName(name)

    for idx, eachCoord in enumerate(coords):
        points.InsertNextPoint(eachCoord[0], eachCoord[1], 0.0)
        point_data.InsertNextTuple1(vals[idx])
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(point_data)

    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(polydata)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputConnection(delaunay.GetOutputPort())
    writer.SetFileName(work_dir)
    writer.Write()

def func(x, t):

    v = 20.0 # um/hour
    t_50 = 16.0 # hours
    rho_max = 0.1
    k = 1.0/20.0
 
    x_0 = v*t
    rho = rho_max*(t/(t_50+t))

    z = rho/(1.0 + np.exp(k*(x-x_0)))
    
    return z
    
def evaluate_fit(coords, v, rho_max, t_50, k):
    
    z_fit = []
    for eachCoord in coords:
        t = eachCoord[1]
        x = eachCoord[0]
        
        x_0 = v*t
        rho = rho_max*(t/(t_50+t))
        z = rho/(1.0 + np.exp(k*(x-x_0)))        
        z_fit.append(z)
    return z_fit

def simple_evaluate_fit(coords, v, rho_max, k, a, beta):
    
    z_fit = []
    for eachCoord in coords:
        t = eachCoord[1]/72.0
        x = eachCoord[0]/1000.0
        
        x_0 = v*t
        x_bar = x-x_0
        rho = a*t*t
        right = 1.0/(1.0 + np.exp(k*(x-x_0)))
        left = beta*t/(1.0 + np.exp(k*(x-x_0+0.1)))
        z = rho*(right-left)
        z_fit.append(z)
    return z_fit

if __name__ == "__main__":
    
    #file_handler = chaste.core.OutputFileHandler("Python/Cornea/ParamSweep_WithConsumption/ParamName_sproutingprobability/ParamValue_0", False)
    file_handler = chaste.core.OutputFileHandler("Python/Cornea/TestSimulationFixedGradient/", False)
    work_dir = file_handler.GetOutputDirectoryFullPath()
    domain_types = ["Planar_2D"]
    
    domain_types = ["Planar_2D", "Planar_3D",  "Circle_2D", "Circle_3D", "Hemisphere"]
    #output_params = ["Line_density", "Tip_density", "Branch_density"]
    #output_params = ["Line_density"]
    output_params = ["Line_density", "Tip_density"]
    
    results_file = open("/home/grogan/results.txt", "w")
    for eachDomainType in domain_types:
        for eachOutputParam in output_params:
            results_dir = work_dir + "/DomainType_" + eachDomainType.replace(" ", "")+"/Run_" + str(0)+"/"
            results_dir += "Sampled_" + eachOutputParam
            x, values  = cornea.postprocessing.sampled_quantity.process_csv(results_dir + ".txt")    
        
            coords = []
            x_eval = []
            t_eval = []
            z_eval = []
            z_simp = []
            for idx in range(2, int(len(x))):
                for jdx in range(0, len(values)):            
                    x_eval.append(x[idx])
                    t = values[jdx][0]
                    t_eval.append(t)
                    coords.append([x[idx], t])
                    results = values[jdx][1]
                    z_simp.append(func(x[idx], t))
                    z = results[idx]
                    z_eval.append(z)
            
            Z, X, T = variables('Z, X, T')
            #Z, X = variables('Z, X')
             
            #v = Parameter(value = 20.0, min=0.0, max=100.0)
            #rho_max = Parameter(value = 0.025, min=0.0, max=0.1)
            #t_50 = Parameter(value = 16.0, min=0.1, max=20.0)
            #k = Parameter(value = 1/20.0, min=1/200.0, max=1.0)
             
            #v, rho_max, t_50, k = parameters('v, rho_max, t_50, k')
            #model = {Z: rho_max*(T/(t_50+T))/(1.0+exp(k*(X-v*T)))}
            #fit = Fit(model, X=np.array(x_eval), T=np.array(t_eval), Z=np.array(z_eval))
             
            v = Parameter(value = 2.0, min=0.0, max=10.0)
            k = Parameter(value = 2.0, min=0.1, max=40.0)
            ap = Parameter(value = 0.1, min=0.01, max=5.0)
            #rho_0 = Parameter(value = 1.0, min=0.0, max=100.0)
            beta = Parameter(value = 1.0, min=-1.0, max=1.0)
            
            #model = {Z: ((rho_0 * T/(ap+T))/(1.0+exp(k*(X-v*T))))}
            model = {Z: (ap*T*T)*(1.0/(1.0+exp(k*(X-v*T))) - (beta*T/(1.0+exp(k*(X-v*t+0.1)))))}

            density_scale_factor = 0.02 # 0.0005
            
            if "2D" in eachDomainType and "Tip" in eachOutputParam:
                density_scale_factor = 0.0005
            elif "3D" in eachDomainType and "Line" in eachOutputParam:    
                density_scale_factor = 0.0001
            elif "3D" in eachDomainType and "Tip" in eachOutputParam:    
                density_scale_factor = 0.0000005 # 0.0005
            x_norm = np.array(x_eval)/1000.0
            z_norm = np.array(z_eval)/density_scale_factor
            
            time_scale = 72.0
            t_norm = np.array(t_eval)/time_scale


            fit = Fit(model, X=x_norm, Z=z_norm, T=t_norm)
            fit_result = fit.execute()
             
            print(fit_result)
            k_res = fit_result.value(k)
            v_res = fit_result.value(v)
            #rho_fit = fit_result.value(rho_0)
            
            rho_fit = 1.0
            a_fit = fit_result.value(ap)
            beta_fit = fit_result.value(beta)
            v_fit = fit_result.value(v)
            #rho_max_fit = fit_result.value(rho_max)
            #t_50_fit = fit_result.value(t_50)
            k_fit = fit_result.value(k)
            
            results_file.write(eachDomainType + ", " + eachOutputParam + ", v: " + str(v_fit) + ", ap: " + str(a_fit) + 
                               ", beta: " + str(a_fit) + ", k: " + str(k_fit) + "\n")
             
            #print fit_result.r_squared
            plot_3d("/home/grogan/old_sampled_" + eachDomainType + "_" + eachOutputParam + ".vtp", coords, z_eval, name="Sampled")

            z_fit = simple_evaluate_fit(coords, v_fit, rho_fit, k_fit, a_fit, beta_fit)
            plot_3d("/home/grogan/old_fit_" + eachDomainType + "_" + eachOutputParam + ".vtp", coords, 
                    np.array(z_fit)*density_scale_factor, name="Fit")
            #plot_3d("/home/grogan/old_simp_" + eachDomainType + "_" + eachOutputParam + ".vtp", name="Fit")
    results_file.close()
    #plt.show()
