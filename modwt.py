import stat
import pywt._wavelet_packets
import json
import pywt.data
import numpy as np
import bioread
from math import log2
import matplotlib.pyplot as plt

class WaveletsTransforms: 
    """Performs MODWT 
    (maximum overlap discrete wavelet transform)


    Notes 
    -----
    Maximum Overlap Discrete Transform  [1]

    References
    ----------
    [1] wavelet methods for Time Series Analysis, Percival & Walden 


    """
    def __init__(self, data, J):
        self.data = data
        self.J = J

    def max_scale(self):
        J = int(log2(len(self.data)) + 1)
        if J % 2 != 0:
            J -= 1
        print(J)
        return J

    def get_wavelet_filters(self, path):
        with open(path, 'r') as file : 
            data = json.load(file)
        self.wavelets = data
    
    def modwt(self, filter):
        h = self.wavelets[filter]['h']
        g = self.wavelets[filter]['g']

        print(h, g)

        N = len(self.data)
        V = [self.data]  # list for scaling coefficients. 
        W = []
        L = len(h)

        for scale in range(1,self.J + 1) : 
            W_scale = []
            V_scale = []
            # N = len(V[-1]) # the previous coef length. 
            for t in range(N):
                k = 0
                k = t
                W_scale_t = h[0] * V[-1][k]
                V_scale_t = g[0] * V[-1][k]
                for n in range(1, L): 
                    k = k - 2 ** (scale-1)
                    if k < 0 : # circular shift
                        k = k + N 
                    W_scale_t += h[n] * V[-1][k] # n 
                    V_scale_t += h[n] * V[-1][k] # n
                W_scale.append(W_scale_t)
                V_scale.append(V_scale_t)
            V.append(V_scale)
            W.append(W_scale)
        return V[-1], W
 
    @staticmethod
    def shift_factor_H(j, L): 
        """
        Computes the factor for the circular shift
        
        """
        L_j = (2**j - 1)*(L - 1) + 1
        v_j = None
        if (L / 2) % 2 == 0: 
            v_j = - int(L_j / 2)
        elif L == 10 or L == 18 : 
            v_j = - int(L_j / 2) + 1
        elif L == 14 : 
            v_j = - int(L_j / 2) - 1
        return v_j

    @staticmethod
    def shift_factor_G(j, L): 
        L_j = (2**j - 1)*(L - 1) + 1
        v_j = None
        if (L / 2) % 2 == 0 : 
            v_j = - int(
                ((L_j - 1) * (L - 2)) / (2* (L - 1))
                )
        elif  L == 10 or L == 18 : 
            v_j = - int(
                ((L_j - 1) * L)/(2 * (L - 1))
            )
        elif L == 14 : 
            v_j = - int(
                ((L_j - 1) * (L - 4)) / (2 * (L - 1))
            )
        return v_j
    
    @staticmethod
    def circular_shift(factor, data):
        return np.roll(data, factor)
        

    def plot_modwt(self, V, W):
        fig, axes = plt.subplots(len(W)+2)
        axes[0].set_title("MODWT decomposition")
        
        axes[0].set_xlim(0, len(self.data) - 1)

        print("taille W", len(W))

        for j, W_j in enumerate(W):            
            ax = axes[-j -2 ]
            factor = self.shift_factor_H(j + 1, 8) # make the lenght of the filter to compute automatically
            print(factor)
            W_j_shifted = self.circular_shift(factor, W_j)
            ax.plot(W_j_shifted, 'g')
            ax.set_ylabel("W%d" % (j+1))
            ax.set_xlim(0, len(self.data) - 1)
        
        ax = axes[0     ]
        factor = self.shift_factor_G(len(W), 8)
        V_shifted = self.circular_shift(factor, V)
        ax.plot(V_shifted, 'm')
        ax.set_ylabel("V%d" % (len(W)))
        ax.set_xlim(0, len(self.data) - 1)

        ax = axes[len(W) + 1]
        ax.set_ylabel('ECG Signal')
        ax.plot(self.data, 'r')
        ax.set_xlim(0, len(self.data) - 1)
        plt.show()

    
    def plot_mra(self, S, D):
        fig, axes = plt.subplots(len(D)+2)
        axes[0].set_title("Multiresolution Analysis")
        
        axes[0].set_xlim(0, len(self.data) - 1)
        for j, D_j in enumerate(D):         
            
            ax = axes[-j -2]
            ax.plot(D_j, 'g')
            ax.set_ylabel("D%d" % (j+1))
            ax.set_xlim(0, len(self.data) - 1)
        
        ax = axes[0]
        ax.plot(S, 'm')
        ax.set_ylabel("S%d" % (len(D)))
        ax.set_xlim(0, len(self.data) - 1)

        ax = axes[len(D) + 1]
        ax.set_ylabel('ECG Signal')
        ax.plot(self.data, 'r')
        ax.set_xlim(0, len(self.data) - 1)
        plt.show()


    
    def imodwt(self, V, W, j, filter):
        """
        For one iteration takes V_J0 and compute V_J0 - 1
        """
        V_j_minus_1 = []
        h = self.wavelets[filter]['h']
        g = self.wavelets[filter]['g']
        N = len(V)
        L = len(h)
        V_j_minus_1
        for t in range(N):
            k = t
            V_j_minus_1_k = h[0] * W[k] + g[0] * V[k]
            for n in range(1, L) :
                k = k + 2 ** (j-1)
                if k >= N : # circular shift
                    k = k - N 
                V_j_minus_1_k += h[n] * W[k] + g[n] * V[k]
            V_j_minus_1.append(V_j_minus_1_k)
        return V_j_minus_1

    def from_V_to_X(self, V_J, W):
        J0 = len(W)
        all_V_j_minus_1 = [
            V_J
        ]
        for scale in range(len(W) -1 , -1, -1):
            print(scale + 1)
            all_V_j_minus_1.append(
                self.imodwt(
                    all_V_j_minus_1[-1],
                    W[scale],
                    scale + 1
                )
            )
        return all_V_j_minus_1[-1]

    
    def mra(self, V, W, J, filter):
        
        D = [] 
        S = []

        # compute Dj
        O_j = [0 for i in range(len(V))]
        # for scale in range(J- 1, -1, -1) : 
        for scale in range(J): 
            W_scale = W[scale]
            
            D_scale = []
            D_scale.append(
                    self.imodwt(
                        O_j , 
                        W_scale , 
                        J, 
                        filter
                    )
                )

            for j in range(scale): 
                D_scale.append(
                    self.imodwt(
                        D_scale[-1] , 
                        O_j , 
                        j + 1, 
                        filter
                    )
                )
            D.append(
                D_scale[-1] #the last result is D_j
            )
        
        # Let's compute S_j. 
        S_scale = [V]
        for scale in range (J - 1, -1, -1):
            S_scale.append(
                self.imodwt(
                    S_scale[-1],
                    O_j, 
                    scale + 1, 
                    filter
                )
            )

        S = S_scale[-1]
        return S, D


if __name__ == '__main__':


    path = "/home/julie/Documents/IRBA/data/Fichier acknoledge centrifugeuse 2009 à 2021/IAM2022/S3 2022/IAM15030_17012022.acq"
    file = bioread.read(path)

    data = file.channels[5].data

    
    ecg = data[:3000]
    print(len(ecg))
    m = WaveletsTransforms(ecg, 4)
    m.get_wavelet_filters("/home/julie/Documents/IRBA/dev/initial_approaches/wavelets/wavelets.json")

    V, W = m.modwt('db6')

    
    m.plot_modwt(V, W)

    S, D = m.mra(V, W, 4, 'db6')

    m.plot_mra(S, D)
    

