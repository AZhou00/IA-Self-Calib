#Interpolation functions - 1D and 2D

def interpolate_1d(y_old,x_old,x_new,method='cubic'):
    """
    y_old is the function value f[x](1D array) evaluated on x=x_old(1D array).
    
    return interpolated f[x] evaluated on x = x_new.
    """
    import numpy as np
    from scipy.interpolate import interp1d
    if (min(x_new) < min(x_old)) or (max(x_new) > max(x_old)):
        raise ValueError('new x axis extends beyond the old x axis')
        
    f = interp1d(x_old, y_old, kind = method)
    y_new = f(x_new)
    
    return y_new

def interpolate_2d(z_old,y_old,x_old,y_new,x_new,method='cubic'):
    """
    z_old is the function value f[y,x](2D array) evaluated on y=y_old(1D array), x=x_old(1D array).
    
    return interpolated f[y,x] evaluated on y=y_new, x=x_new.
    """
    import numpy as np
    from scipy.interpolate import griddata
    if (min(y_new) < min(y_old)) or (max(y_new) > max(y_old)):
        raise ValueError('new y axis extends beyond the old y axis')
    if (min(x_new) < min(x_old)) or (max(x_new) > max(x_old)):
        raise ValueError('new x axis extends beyond the old x axis')
        
    #converts the original z.shape=(leny,lenx) format to z as a list of points (y0,x0),(y1,x1)
    mesh = np.array(np.meshgrid(y_old, x_old))
    points = mesh.T.reshape(-1, 2)
    values = z_old.flatten()
    
    y_new_mesh, x_new_mesh = np.meshgrid(x_new,y_new)
    z_new = griddata(points, values, (x_new_mesh,y_new_mesh), method=method)
    
    return z_new
