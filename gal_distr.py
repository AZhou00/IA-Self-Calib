def get_gal_distr(z_bin_num,deltaz,z_res,z_max,model,outputpath,REMOTE_COSMOSIS_FLAG):
    
    """
    model: string. e.g. "z2exp_gausPDF"
    z_max = z_bin_num*deltaz+constant(~2)
    
    bin1 from       0 - 1deltaz
    bin2 from 1deltaz - 2deltaz
    bin3 from 2deltaz - 3deltaz
    ...
    
    generate galaxy distributions according to photo_n_distr convolved with photo_z_PDF 
        function values all evaluated at z_mid for each redshift histogram bins.
    
    generate BINNED photo-z AND true-z distributions with data structure 
        dtype=[('Z_LOW', '>f8'), ('Z_MID', '>f8'), ('Z_HIGH', '>f8'), ('BIN1', '>f8'), ('BIN2', '>f8'), ...]
    
    genearte AGGREGATE (from z_min to z_max) photo-z AND true-z distributions
    
    
    """
    import os,sys
    import numpy as np
    import scipy as sp
    import scipy.integrate as intg
    import io_fncs
    from matplotlib import pyplot as plt
    
    if REMOTE_COSMOSIS_FLAG == True:
        import fitsio as fio
    
    hist_bin_num = int(z_max*z_res) #how many histogram bins there are as opposed to number of redshift bins    
    
    data_type = [('Z_LOW', '>f8'), ('Z_MID', '>f8'), ('Z_HIGH', '>f8')]
    for i in range(z_bin_num):
        data_type.append(('BIN%i'%(i+1), '>f8'))
        
    photo_z_file = np.empty(hist_bin_num,dtype=data_type)
    true_z_file  = np.empty(hist_bin_num,dtype=data_type)
    for i in range(hist_bin_num):
        z_low  = i*1/z_res
        z_high = z_low + 1/z_res
        z_mid  = (z_low+z_high)/2
        photo_z_file[i][0] = z_low
        photo_z_file[i][1] = z_mid
        photo_z_file[i][2] = z_high 
        true_z_file[i][0]  = z_low
        true_z_file[i][1]  = z_mid
        true_z_file[i][2]  = z_high 
    
    if model == "z2exp_gausPDF":
        gaussian_spread = 0.05
        
        def photo_n_distr(z):
            return z*z*np.exp(-z/0.5) 

        z_list = photo_z_file['Z_MID']
        n_p_list = photo_n_distr(z_list)
        for i in range(z_bin_num):
            
            #photo-z
            bin_z_low = i*deltaz
            bin_z_high = bin_z_low+deltaz
            selection = 1*((z_list-1/(4*z_res))>=bin_z_low)*((z_list+1/(4*z_res))<=bin_z_high) 
            #use 1/4 to account for rounding errors
            #can check by print out selection. will see true's for bins are disjoint and complete.
            #print(selection)
            n_p_bin_list = n_p_list*selection
            n_p_bin_list = n_p_bin_list/(np.sum(n_p_bin_list)*1/z_res) #normalization
            photo_z_file['BIN%i'%(i+1)] = n_p_bin_list
            
            #true-z
            n_true_list = np.zeros(hist_bin_num)
            for j,val in enumerate(n_p_bin_list):
                if val != 0: #if val=0, then because of the sharp edge for the n_p_bin distribution, it does not contribute to n(true z)
                    zphoto = z_list[j]
                    std = gaussian_spread*(1+zphoto)
                    prefector = 1/(np.sqrt(2*np.pi)*std)
                    n_true_list += n_p_bin_list[j]*prefector*np.exp(-0.5*((z_list-zphoto)/std)**2)
            n_true_list = n_true_list/(np.sum(n_true_list)*1/z_res)
            true_z_file['BIN%i'%(i+1)] = n_true_list
    
    #IO
    io_fncs.create_dir(os.path.join(outputpath,'n_z'))
    
    io_fncs.save_config_file(os.path.join(outputpath,'n_z/values.txt'),
                 ['deltaz',deltaz,
                  'z_bin_num',z_bin_num,
                  'z_max',z_max,
                  'z_res',z_res,
                  'z_res_spacing',1/z_res,
                  'model',model,
                  'dtype',str(data_type)])

    
    txt_format = '%f,%f,%f'
    for i in range(z_bin_num):
        txt_format += ',%f'
    np.savetxt(os.path.join(outputpath,'n_z/n_photo_z.txt'), photo_z_file, fmt=txt_format)
    np.savetxt(os.path.join(outputpath,'n_z/n_true_z.txt'), true_z_file, fmt=txt_format)
    if REMOTE_COSMOSIS_FLAG == True:
        #writing fits file
        #clobber=True ->overwrite
        #extname -> set hduname
        #fits source code at
        #https://github.com/esheldon/fitsio/blob/master/fitsio/fitslib.py
        fio.write(os.path.join(outputpath,'n_z/n_photo_z.fits'),photo_z_file,clobber=True,extname='NZ_SAMPLE')
        fio.write(os.path.join(outputpath,'n_z/n_true_z.fits'),true_z_file,clobber=True,extname='NZ_SAMPLE')
        
    plt.subplots(figsize=(8,6))
    for i in range(z_bin_num):
        plt.plot(z_list, photo_z_file['BIN%i'%(i+1)],label = 'PHOTOZ BIN%i'%(i+1))
        plt.plot(z_list, true_z_file['BIN%i'%(i+1)],linestyle='dashed',label = 'TRUEZ  BIN%i'%(i+1))
    plt.xlabel('redshift z')
    plt.ylabel('galaxy redshift distribution')
    plt.title('n(z), normalized for each bin seperately')
    #plt.legend(bbox_to_anchor=(1, 1.05))
    plt.tight_layout()
    plt.savefig(os.path.join(outputpath,'n_z/combined_n_z.png'))
    #plt.show()
    plt.close()
    
    return [photo_z_file,true_z_file]
