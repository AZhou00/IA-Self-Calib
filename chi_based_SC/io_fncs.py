
def create_dir(new_folder_path):
    """
    create a directory if it does not exist.
    """
    import os
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
    return

def append_to_file(filename,row,delimeter=' '):
    """
    append row to the file 'filename'. 
    'row' should be a numpy array.
    no end-of-the-line delimeter
    
    """
    import numpy as np
    with open(filename,'a') as f:
        for index,el in enumerate(row):
            if index != (len(row)-1):
                f.write(str(el)+delimeter)
            else:
                f.write(str(el))
            #no end of the line delimeter
        f.write("\n")
    return

def load_num_file(filename,delimeter=' '):
    """
    load a numeric txt file. floats only.
    """
    with open(filename,'rb') as f:
        out = np.loadtxt(f,delimiter = delimeter) 
    return out

def save_config_file(filename,var_names):
    """
    var_names: ['var1',var1,'var2',var2,'var3',var3,...]
    write to files:
    'var1' = var1
    'var2' = var2
    ...
    - delete pre-existing file if any.
    """
    import numpy as np
    with open(filename,'w') as f:
        if len(var_names)%2 != 0: raise Exception('variable names not paired correctly')
        index = 0
        for index in range(len(var_names)):
            if index%2 == 0:
                name  = str(var_names[index])
                value = str(var_names[index+1])
                f.write(name+' = '+value)
                f.write("\n")
    return
