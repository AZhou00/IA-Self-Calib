def get_basic_cosmo_lensing(d_horiz,Omega0,OmegaLambda,Omegak,z_list,outputpath,USECOSMOSIS = False,data_path = ''):
    """
    if USECOSMOSIS is on (by default.) then chi is read from cosmosis output, and then interpolated on z_list.
    
    save chi(the cmv distance in c/H0 unit)
    save W_L(the lensing kernel, with chi in the above mentioned unit)
    save their plots
    
    functions of z_list, mid bin, bins are as defined in n_z
    
    """
    import os,sys
    import numpy as np
    import scipy as sp
    import scipy.integrate as intg
    import io_fncs
    from matplotlib import pyplot as plt
    from scipy import interpolate
    
    if (USECOSMOSIS == True and data_path == ''):
        raise Exception('if using chi from cosmosis, need to specify the input folder path. The input path folder should both contain chi and the corresponding z_list to use.')
    
    io_fncs.create_dir(os.path.join(outputpath,'basic_fncs'))
    if USECOSMOSIS == False:
        #comoving distance function, in unit of Mpc (since d_horiz in unit of Mpc)  ###(CHECKED)###
        #Troxel and Ishak 2014 equation (8) without c=1.
        chi_list = np.array([],dtype='float64')
        def chi_integrand(z):
            return 1/np.sqrt(Omega0*(1+z)**3+Omegak*(1+z)**2+OmegaLambda)
        def chi(z):
            return intg.quad(chi_integrand,0,z)[0]
        for z in z_list: chi_list = np.append(chi_list,chi(z))
        chi_list = chi_list*d_horiz #in Mpc
        plt.plot(z_list, chi_list,color='red')
        
    if USECOSMOSIS == True:
        with open(os.path.join(data_path,'d_m.txt'),'rb') as f:
            chi_cosmosis = np.loadtxt(f) #already in unit of Mpc by default
        with open(os.path.join(data_path,'z.txt'),'rb') as f:
            z_cosmosis = np.loadtxt(f)
        chi_func = interpolate.interp1d(z_cosmosis,chi_cosmosis,kind='linear') 
        chi_list = chi_func(z_list)
        plt.plot(z_list, chi_list,color='black')
        
    
    plt.xlabel('redshift z')
    plt.ylabel('comoving distance [Mpc] chi(z) \n c/H0= %1.2f, Omega0=%1.2f, OmegaLambda=%1.2f' %(d_horiz, Omega0,OmegaLambda))
    plt.tight_layout()
    plt.savefig(os.path.join(outputpath,'basic_fncs/chi.png'))
    plt.close()
    
    
    # Lensing Kernel List  ###(CHECKED)###
    # Troxel Ishak eqn 38 https://arxiv.org/pdf/1407.6990.pdf
    # W_L_list[z_G,z_L], shape(len(z_list),len(z_list))
    W_L_list = np.empty((len(z_list),len(z_list)),dtype='float64')
    for z_G_index in range(len(z_list)):
        chi_G = chi_list[z_G_index]
        W_L_list[z_G_index] = (1/d_horiz**2)*1.5*Omega0*(1+z_list)*chi_list*(1-chi_list/chi_G)
    W_L_list[W_L_list<0]=0 #set negative entries to 0 (eqv. set chi_L>chi_G entries to 0)
    #test = W_L_list.clip(min=0) This is slightly slower than the above method
    
    test_redshift_index = int(1.2/(z_list[1]-z_list[0]))
    plt.plot(z_list, W_L_list[test_redshift_index])
    plt.xlabel('redshift z')
    plt.ylabel('lensing kernel W_L [Mpc] with c/H0= %1.2f and source z_G=%1.4f ' %(d_horiz, z_list[test_redshift_index]))
    plt.tight_layout()
    plt.savefig(os.path.join(outputpath,'basic_fncs/W_L.png'))
    plt.close()
    
    return [chi_list,W_L_list]

def get_integrated_lensing_kernel(W_L_list,true_z_file,z_list,z_bin_num,outputpath):  ###(CHECKED)###
    """
    get integrated lensing kernel W_i for each redshift bin i, i = 1,2,3,4...
    
    Flat geometry
    
    return W_i[z_bin_#, z_L], size = (z_bin_num,len(z_list))
    
    assume uniform z-intervals
    z_list should = true_z_file['Z_MID'] = photo_z_file['Z_MID'] for consistency
    z_bin_num should = len(true_z_file[0])-3 = len(photo_z_file[0])-3 for consistency 
        (3 comes from the first 3 columns Z_LOW,Z_MID,Z_HIGH for each true_z_file entry)
    """
    import os,sys
    import numpy as np
    import io_fncs
    from matplotlib import pyplot as plt
    
    if not np.array_equal(z_list,true_z_file['Z_MID']): raise Exception('z_list failed consistency test')
    if len(W_L_list[0]) != len(z_list): raise Exception('W_L_list size not consistent with imput z_list')
    if z_bin_num  != len(true_z_file[0])-3: raise Exception('z_bin_num failed consistency test')
    if abs((z_list[1]-z_list[0])/(z_list[-1]-z_list[-2])-1) >= 0.01: raise Exception('nonuniform z_list detected >=1%')
        
    io_fncs.create_dir(os.path.join(outputpath,'basic_fncs'))
    
    z_spacing = z_list[1]-z_list[0]
    W_i = np.empty((z_bin_num,len(z_list)))
    
    plt.subplots(figsize=(6,6))
    for bin_index in range(z_bin_num):
        n_bin = true_z_file['BIN%i'%(bin_index+1)]
        for z_L_index in range(len(z_list)):
            #flat geometry
            W_i[bin_index,z_L_index] = sum(W_L_list[:,z_L_index]*n_bin)*z_spacing#sum((W_L_list.T)[z_L_index]*n_bin)*z_spacing #using the true redshift distribution
        plt.plot(z_list,W_i[bin_index],label='bin_i=%i'%(bin_index+1))

    plt.xlabel('redshift z')
    plt.ylabel('integrated kernel W_i')
    plt.legend(bbox_to_anchor=(1, 0),loc='lower left') #bbox(x,y)
    plt.tight_layout()
    plt.savefig(os.path.join(outputpath,'basic_fncs/W_i.png'))
    #plt.show()
    plt.close()
    return W_i
