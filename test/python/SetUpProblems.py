
import chaste
import microvessel_chaste
import microvessel_chaste.pde
from microvessel_chaste.utility import *

def set_up_reference_scales():
    
    reference_concentration = 1.e-9 *mole_per_metre_cubed()
    reference_length_scale = 1.e-6 * metre()
    reference_time_scale = 3600.0 * second()
    BaseUnits.Instace().SetReferenceLengthScale(reference_length_scale)
    BaseUnits.Instace().SetReferenceTimeScale(reference_time_scale)
    BaseUnits.Instace().SetReferenceConcentrationScale(reference_concentration)
    
def get_pde(dimension=2):
    
    if dimension ==2:
        pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde2()
    else:
        pde = microvessel_chaste.pde.CoupledVegfPelletDiffusionReactionPde3()
    
    
    
        boost::shared_ptr<CoupledVegfPelletDiffusionReactionPde<DIM> > p_pde = CoupledVegfPelletDiffusionReactionPde<DIM>::Create();
    units::quantity<unit::diffusivity> vegf_diffusivity(6.94e-11 * unit::metre_squared_per_second);
    units::quantity<unit::rate> vegf_decay_rate((-0.8/3600.0) * unit::per_second);
    p_pde->SetIsotropicDiffusionConstant(vegf_diffusivity);
    p_pde->SetContinuumLinearInUTerm(vegf_decay_rate);
    units::quantity<unit::concentration> initial_vegf_concentration(3.93e-1*unit::mole_per_metre_cubed);
    p_pde->SetCurrentVegfInPellet(initial_vegf_concentration);
    p_pde->SetPelletBindingConstant(100.0);
    p_pde->SetPelletDepth(50.0e-6*unit::metres);
    p_pde->SetPelletSurfaceArea(2000e-6*unit::metres*p_pde->GetPelletDepth());
    p_pde->SetCorneaPelletPermeability(0.002*p_pde->GetCorneaPelletPermeability());
    return p_pde;

