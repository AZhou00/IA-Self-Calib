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
    import numpy as np
    with open(filename,'rb') as f:
        out = np.loadtxt(f,delimiter = delimeter) 
    return out
