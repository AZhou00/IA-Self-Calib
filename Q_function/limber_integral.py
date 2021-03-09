def limber_integral(Del_mg,limber_weight,chi_list,k_list,z_list,ell_list):
    """
    - compute Integral(Del_mg(k=l/chi(z))*limber_weight)dz
    - all functions of z takes value on z_list (need to interpolate first)
    - z_list must be uniform
    - this method uses interpolation (default:cubic)
    
    """
    import numpy as np
    from scipy import interpolate

    if abs((z_list[1]-z_list[0])/(z_list[-1]-z_list[-2])-1) >= 0.01: raise Exception('nonuniform z_list detected >=1%')
    
    z_interval = z_list[1]-z_list[0]
    C_l = np.zeros(len(ell_list))
    
    for z_index in range(len(z_list)):
        temp_integrand = z_interval*Del_mg[z_index]*limber_weight[z_index]
        temp_ell_coord = k_list*chi_list[z_index]
        C_l_increment_function = interpolate.interp1d(temp_ell_coord,temp_integrand,bounds_error=False,fill_value=0,kind='linear')
        C_l += C_l_increment_function(ell_list)

    return 2*np.pi*np.pi*C_l/(ell_list**3)
