import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageColor, ImagePalette
import os
from datetime import datetime
import glob

def radiance_to_bt(RAD,fk1,fk2,bc1,bc2):
    T=fk2/(np.log((fk1/RAD)+1))
    ir_BT=bc2*T+bc1
    return ir_BT

def radiance_to_ref(RAD,d,ensu):
    return (RAD * np.pi * d) /ensu 

def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    m = np.arange(0,100,1).reshape((10,10))
    n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,
                                                   ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
    return ndarray

def get_like_index(array,search_val):
    aux= abs(array-search_val).tolist()
    mindist= min(aux)
    indice=aux.index(mindist)
    return indice

def creapaleta(i_tempcol):
    color = open("colorHEX.txt","r")
    cmap = color.readlines()
    lut = []
    for i in range(256):
        # If the index is above the minimum color, then we change the color to the one in the $
        if (i > i_tempcol):
            # WHY MULTIPLIED BY 360
            sc = int((i - i_tempcol)*255/(255 - i_tempcol))
            color = ImageColor.getrgb(cmap[sc])
        
        else:
    # If the index is below we simply leave the gray index
            c = int(255.0*i/i_tempcol)
            color = [c,c,c]
        lut.extend(color)
    return lut

def temptoindex(t, temp_max, temp_min):
    return (t - temp_min)*255/(temp_max - temp_min + .0001)

def get_wind_file(date_string,path):
    # Convert the string to a datetime object
    date_object = datetime.strptime(date_string, "%Y.%m%d.%H%M")

    # Calculate the start of the year
    start_of_year = datetime(date_object.year, 1, 1)

    # Calculate the difference in days
    day_of_year = (date_object - start_of_year).days + 1
    print("%%%%%%%day_of_year :  ",day_of_year)
    pattern=path+"OR_ABI-L2-DMWF-M6C02_G16_s"+ date_object.strftime("%Y")+str(day_of_year).zfill(3)+date_object.strftime("%H")+"*"
    print("%%%%%%%Patte :  ",pattern)
    files = glob.glob(pattern)

    # Print the files that match the pattern
    for file in files:
      return file
    return ""
