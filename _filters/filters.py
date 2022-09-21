from audioop import reverse
from math import sqrt
import os
import json

class Filter: 
    """
    Builds filter for MODWT. 

    Parameters
    ----------

    self.name : str.
        name of the filter
    self.g : list.
        scaling coefficients
    self.h : list.
        wavelet coefficients
    
    """
    def __init__(self, filter_name):
        """
        Parameters
        ----------
        filter_name : str.
            name of filter. 
            Need to be choosen among : haar, db4, db6,
        """

        self.name = filter_name
        with open('_filters/wavelets.json', 'r') as f :
            self.dict = json.load(f)
        try : 
            self.g = [coef / sqrt(2)
                for coef in self.dict[filter_name]['g'] 
            ]
        
            self.L = len(self.g)
            self.h = Filter.qmf(self.g, self.L)
        except : 
            print(filter_name + " not in the filter bank")
            raise

    @staticmethod
    def qmf(g, L): 
        """
        Quadrature miror filter.
        Computes h from g. 

        Parameters
        ----------
        g : list.
            scaling coefficients.
        L : int.
            size of the coefficients. 

        Returns
        -------
        h : list.
            wavelet coefficients.
        """
        h = [ 
            g[l] * (-1) ** (l - 1)
            for l in range(L - 1, -1, -1)
        ]
        return h

    def __str__(self):
        res = self.name
        g = ','.join(
            [
                str(c) for c in self.g
            ]
        )
        h = ','.join(
            [
                str(c) for c in self.h
            ]
        )
        res += "\n Wavelet Coefficients " + h
        res += "\n Scaling Coefficients " + g
        return res
    