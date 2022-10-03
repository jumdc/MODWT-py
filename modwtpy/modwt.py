
import numpy as np
import matplotlib.pyplot as plt


class MODWT: 
    def __init__(self, data, filter, J): 
        """
        Paramters
        ---------
        data : list.
            data to perform the MODWT on.
        filter : Filter.
            filter to use during the MODWT. 
        J : int.
            level of the MODWT 
            i.e : number of scales.
        """
        self.data = data
        self.filter =  filter
        self.L = len(filter.g)
        self.J = J
        self.N = len(data)
        self.W = []
        self.V = []

    def modwt(self): 
        """ Computes MODWT. 
        Notes 
        -----
        Maximum Overlap Discrete Transform  [1]

        Returns
        -------
        Wavelet and scaling coefficients. 

        References
        ----------
        [1] wavelet methods for Time Series Analysis, Percival & Walden 
        """
        h = self.filter.h
        g = self.filter.g

        N = len(self.data)
        V = [self.data]  # list for scaling coefficients. 
        W = []

        for scale in range(1,self.J + 1) : 
            W_scale = []
            V_scale = []
            # N = len(V[-1]) # the previous coef length. 
            for t in range(N):
                k = t
                W_scale_t = h[0] * V[-1][k]
                V_scale_t = g[0] * V[-1][k]
                for n in range(1, self.L): 
                    k = k - 2 ** (scale-1)
                    if k < 0 : # circular shift
                        k = k + N 
                    W_scale_t += h[n] * V[-1][k] # n 
                    V_scale_t += g[n] * V[-1][k] # n
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
        """
        Plots MODWT wavelet and scaling coefficients. 
        """
        fig, axes = plt.subplots(len(W)+2)
        axes[0].set_title("MODWT decomposition")
        
        axes[0].set_xlim(0, len(self.data) - 1)


        for j, W_j in enumerate(W):            
            ax = axes[-j -2 ]
            factor = self.shift_factor_H(j + 1, self.L) # make the lenght of the filter to compute automatically
            W_j_shifted = self.circular_shift(factor, W_j)
            ax.plot(W_j_shifted, 'g')
            ax.set_ylabel("W%d" % (j+1))
            ax.set_xlim(0, len(self.data) - 1)
        
        ax = axes[0]
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

    def modwt_backward(self, V, W, j):
        """
        For one iteration takes V_J0 and compute V_J0 - 1
        """
        V_j_minus_1 = []
        h = self.filter.h
        g = self.filter.g
        for t in range(self.N):
            k = t
            V_j_minus_1_k = h[0] * W[k] + g[0] * V[k]
            for n in range(1, self.L) :
                k = k + 2 ** (j-1)
                if k >= self.N : # circular shift
                    k = k - self.N 
                V_j_minus_1_k += h[n] * W[k] + g[n] * V[k]
            V_j_minus_1.append(V_j_minus_1_k)
        return V_j_minus_1

    def imodwt(self, V, W):
        all_V_j_minus_1 = [
            V
        ]
        for scale in range(self.J -1 , -1, -1):
            all_V_j_minus_1.append(
                self.modwt_backward(
                    all_V_j_minus_1[-1],
                    W[scale],
                    scale + 1
                )
            )
        return all_V_j_minus_1[-1]

    def mra(self, V, W):
        D = [] 
        S = []

        # compute Dj
        O_j = [0 for i in range(len(V))]
        # for scale in range(J- 1, -1, -1) : 
        for scale in range(self.J): 
            W_scale = W[scale]
            
            D_scale = []
            D_scale.append(
                    self.modwt_backward(
                        O_j , 
                        W_scale , 
                        self.J
                    )
                )
            for j in range(1, scale): 
                D_scale.append(
                    self.modwt_backward(
                        D_scale[-1] , 
                        O_j , 
                        j + 1
                    )
                )
            D.append(
                D_scale[-1] #the last result is D_j
            )
        
        # Let's compute S_j. 
        S_scale = [V]
        for scale in range (self.J - 1, -1, -1):
            S_scale.append(
                self.modwt_backward(
                    S_scale[-1],
                    O_j, 
                    scale + 1
                )
            )
        S = S_scale[-1]
        return S, D

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
