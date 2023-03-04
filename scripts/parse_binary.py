import numpy as np
from struct import *
import zlib

def unpack_array(bin_array, numpeaks, array_type):

    if array_type == 'peaks':
        parsed_array = np.array([], dtype=float)
    elif array_type == 'intensity': 
        parsed_array = np.array([], dtype='f')

    if array_type == 'peaks':
        a_format = "<" + numpeaks + 'd'
    elif array_type == 'intensity':
        a_format = "<" + numpeaks + 'f'
    
    try:
        decompress_bin_array = zlib.decompress(bin_array)
        decoded_array = unpack(a_format, decompress_bin_array)
        parsed_array = np.asarray(decoded_array)
    except:
        try:
            decoded_array = unpack(a_format, bin_array)
        except:
            if "d" in a_format:
                new_format = "<" + str(len(bin_array)) + "d"
            elif "f" in a_format:
                new_format = "<" + str(len(bin_array)) + "f"
                                       
            decoded_array = unpack(new_format, bin_array)
            
        parsed_array = np.asarray(decoded_array)
    return parsed_array
