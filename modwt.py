import stat
import pywt._wavelet_packets
import pywt.data
import numpy as np
import bioread
from math import log2
import matplotlib.pyplot as plt

class MODWT: 
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
    
    def modwt(self):
        h = [
            -0.0105974018, 0.0328830117, 0.0308413818, 
            -0.1870348117, -0.0279837694, 0.6308807679, 
            0.7148465706, 0.2303778133
            ]
        g = [
            -0.2303778133, 0.7148465706, -0.6308807679,
            -0.0279837694, 0.1870348117, 0.0308413818, 
            -0.0328830117, -0.0105974018
        ]

        N = len(self.data)
        V = []  # list for scaling coefficients. 
        V.append(self.data)
        W = []
        L = len(h)

        for scale in range(1,self.J + 1) : 
            print("scale", scale)
            W_scale = []
            V_scale = []
            # N = len(V[-1]) # the previous coef length. 
            for t in range(N):
                k = 0
                k = t
                W_scale_t = h[0]*V[-1][k]
                V_scale_t = g[0]*V[-1][k]
                for n in range(1, L ): 
                    k = k - 2**(scale-1)
                    if k < 0 : # circular shift
                        k = k + N 
                    W_scale_t += h[n]*V[-1][k] # n 
                    V_scale_t += h[n]*V[-1][k] # n
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
        print(data[:10], data[-10:])
        res = np.roll(data, factor)
        print(res[:10], res[-10:])
        return np.roll(data, factor)
        

    def plot_modwt(self, V, W):
        fig, axes = plt.subplots(len(W)+2)
        axes[0].set_title("MODWT decomposition")
        
        axes[0].set_xlim(0, len(self.data) - 1)

        print("taille W", len(W))

        for j, W_j in enumerate(W):            
            ax = axes[j]
            factor = self.shift_factor_H(j + 1, 8) # make the lenght of the filter to compute automatically
            print(factor)
            W_j_shifted = self.circular_shift(factor, W_j)
            ax.plot(W_j_shifted, 'g')
            ax.set_ylabel("W%d" % (j+1))
            ax.set_xlim(0, len(self.data) - 1)
        
        ax = axes[len(W)]
        factor = self.shift_factor_G(len(W), 8)
        V_shifted = self.circular_shift(factor, V)
        ax.plot(V_shifted, 'm')
        ax.set_ylabel("V%d" % (len(W)))
        ax.set_xlim(0, len(self.data) - 1)

        ax = axes[len(W) + 1]
        ax.set_ylabel('ECG Signal')
        ax.plot(self.data, 'r')

        plt.show()


    def mra(self):
        pass

if __name__ == '__main__':


    path = "/home/julie/Documents/IRBA/data/Fichier acknoledge centrifugeuse 2009 à 2021/IAM2022/S3 2022/IAM15030_17012022.acq"
    file = bioread.read(path)

    data = file.channels[5].data

    ecg = data[:1000]
    print(len(ecg))
    m = MODWT(ecg, 6)


    V, W = m.modwt()
    m.plot_modwt(V, W)

 

    # dwt = MODWT(ecg)
    # dwt.max_j()
    # dwt.modwt()
    # ecg = pywt.data.ecg()

    # # set trim_approx to avoid keeping approximation coefficients for all levels

    # # set norm=True to rescale the wavelets so that the transform partitions the
    # # variance of the input signal among the various coefficient arrays.

    # coeffs = pywt.swt(ecg, wavelet='db10', trim_approx=True, norm=True, level = 6)

    # # coeffs = pywt.mra(ecg, wavelet='db4', transform='swt')

    # ca = coeffs[0]
    # details = coeffs[1:]

    # print("Variance of the ecg signal = {}".format(np.var(ecg, ddof=1)))

    # variances = [np.var(c, ddof=1) for c in coeffs]
    # detail_variances = variances[1:]
    # print("Sum of variance across all SWT coefficients = {}".format(
    #     np.sum(variances)))

    # # Create a plot using the same y axis limits for all coefficient arrays to
    # # illustrate the preservation of amplitude scale across levels when norm=True.
    # ylim = [ecg.min(), ecg.max()]

    # fig, axes = plt.subplots(len(coeffs) + 1)
    # axes[0].set_title("normalized SWT decomposition")
    # axes[0].plot(ecg)
    # axes[0].set_ylabel('ECG Signal')
    # axes[0].set_xlim(0, len(ecg) - 1)
    # axes[0].set_ylim(ylim[0], ylim[1])

    # for i, x in enumerate(coeffs):
    #     ax = axes[-i - 1]
    #     ax.plot(coeffs[i], 'g')
    #     if i == 0:
    #         ax.set_ylabel("A%d" % (len(coeffs) - 1))
    #     else:
    #         ax.set_ylabel("D%d" % (len(coeffs) - i))
    #     # Scale axes
    #     ax.set_xlim(0, len(ecg) - 1)
    #     ax.set_ylim(ylim[0], ylim[1])


    # # reorder from first to last level of coefficients
    # level = np.arange(1, len(detail_variances) + 1)

    # # create a plot of the variance as a function of level
    # # plt.figure(figsize=(8, 6))
    # # fontdict = dict(fontsize=16, fontweight='bold')
    # # plt.plot(level, detail_variances[::-1], 'k.')
    # # plt.xlabel("Decomposition level", fontdict=fontdict)
    # # plt.ylabel("Variance", fontdict=fontdict)
    # # plt.title("Variances of detail coefficients", fontdict=fontdict)
    plt.show()