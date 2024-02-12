from cosmosis.datablock import names, option_section as opt
import numpy as np
import matplotlib.pyplot as plt
import camb
import warnings
from scipy.interpolate import interp1d
    
def setup(options):
    return []
    
def execute(block, config):
    # Getting Plin at z=zstar
    Plin_star=block[names.matter_power_lin,"P_K"][-1,:]
   
    # Getting sigma8 at z=zstar
    sig8_star=block[names.growth_parameters,"sigma_8"][-1]
    
    # Computing the effective power spectra that goes into the linear part
    Peff=np.vstack((Plin_star/sig8_star**2,block[names.matter_power_lin,"P_K"][1:,:]))

    # Replacing Plin(k,z=0) by the effective spectra
    block[names.matter_power_lin,"P_K"] = Peff

    # Getting the boost from the datablock
    B=block[names.matter_power_nl,"B"]
    
    # Replacing the effective nl spectra
    Peff2 = B*Plin_star/sig8_star**2
    block[names.matter_power_nl,"P_K"] = Peff2
    
    #Printing bhat fiducials using b=1
    z=block[names.matter_power_nl,"z"]
    print('zstar',z[-1])
    sig8_z=block[names.growth_parameters,"sigma_8"]
    bhat_z=1*sig8_z
    bhat_func=interp1d(z,bhat_z)
    #print('bhat fiducials',bhat_func([0.2946683204029144, 0.46670347176087024, 0.626278460689437, 0.771314100934901,0.9,1.0]))
    print('bhat fiducials:',bhat_func([0.2796,0.4367,0.5795,0.7271,0.8486]))
    
    return 0
