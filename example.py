import sys
from modwtpy.filters import Filter
from modwtpy.modwt import MODWT
# import matplotlib.pyplot as plt


if __name__ == '__main__':


    
    # Perform MODWT 

    # Initiate the filter

    la8 = Filter('la8')
    c2 = [0.2, -0.4, -0.6, -0.5, -0.8, -0.4, -0.9, -0.2, 0.1, -0.1, 0.1, -0.7, 0.9, 0, 0.3]

    
    m = MODWT( # initite MODWT class
        c2, 
        la8, 
        1
    )

    V, W = m.modwt()
    m.plot_modwt(V, W) # plot the Wavelet and Sclaing coefficients

    
    X = m.imodwt(V,W) # Inverse MODWT 


    S, D = m.mra(V, W) # Multiresolution analysis
    m.plot_mra(S, D)
