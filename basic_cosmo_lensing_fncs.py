def get_basic_cosmo_lensing(Omega0,OmegaLambda,Omegak,z_list,outputpath):
    
    import os,sys
    import numpy as np
    import scipy as sp
    import scipy.integrate as intg
    import io_fncs
    from matplotlib import pyplot as plt
    
    """
    save chi(the cmv distance in c/H0 unit)
    save W_L(the lensing kernel, with chi in the above mentioned unit)
    save their plots
    
    functions of z_list, mid bin, bins are as defined in n_z
    
    """
    io_fncs.create_dir(os.path.join(outputpath,'basic_fncs'))
    
    #comoving distance function, in unit of c/H0
    #Troxel and Ishak 2014 equation (8) without the factor of H0
    chi_list = np.array([],dtype='float64')
    def chi_integrand(z):
        return 1/np.sqrt(Omega0*(1+z)**3+Omegak*(1+z)**2+OmegaLambda)
    def chi(z):
        return intg.quad(chi_integrand,0,z)[0]
    for z in z_list: chi_list = np.append(chi_list,chi(z))
    plt.plot(z_list, chi_list)
    plt.xlabel('redshift z')
    plt.ylabel('comoving distance chi(z) \n Omega0=%1.2f, OmegaLambda=%1.2f' %(Omega0,OmegaLambda))
    plt.savefig(os.path.join(outputpath,'basic_fncs/chi.png'))
    plt.close()
    # Lensing Kernel List
    # W_L_list[z_list,z_list]
    W_L_list = np.empty((len(z_list),len(z_list)),dtype='float64')
    for z_G_index in range(len(z_list)):
        chi_G = chi_list[z_G_index]
        W_L_list[z_G_index] = 1.5*Omega0*(1+z_list)*chi_list*(1-chi_list/chi_G)
    W_L_list[W_L_list<0]=0 #set negative entries to 0 (eqv. set chi_L>chi_G entries to 0)
    #test = W_L_list.clip(min=0) This is slightly slower than the above method
    test_redshift_index = int(1.2/(z_list[1]-z_list[0]))
    plt.plot(z_list, W_L_list[test_redshift_index])
    plt.xlabel('redshift z')
    plt.ylabel('lensing kernel W_L with source z_G=%1.4f' %(z_list[test_redshift_index]))
    plt.savefig(os.path.join(outputpath,'basic_fncs/W_L.png'))
    plt.close()
    return [chi_list,W_L_list]
