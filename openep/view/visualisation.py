'''
Library for importing, exporting and doing visualisation on the triangular meshes
'''

import scipy.io as sio



file_path = '../openep/openep-examples-main/openep_dataset_1.mat'



# Function to read the file
# return the userdata
def load_mat(file_path):
    '''
    Function to read MAT file

    Parameters
    -----------
    file_path : str
        absolute or relative path to the MAT file

    Returns
    ------------
    main_file : dict
        dictionary with variable names as keys and loaded matrices as values
    '''

    try:
        main_file = sio.loadmat(file_path)
        return main_file
    except:
        print('Error: Unable to read file. Check MAT file path')
        return 0

coredata = load_mat(file_path)
