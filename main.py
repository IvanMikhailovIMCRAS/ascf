from ascf import FlatBrush, TwoBrush
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')

if __name__ == "__main__":
    # for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
    #     fb = FlatBrush(N=100, sigma=0.1, chi=chi)
    #     plt.plot(fb.z, fb.phi, label=f"$\chi={chi}$")
    # plt.legend()
    # plt.show()
    
    N = 500
    sigma = 0.1
    colors = {0.0: 'blue', 0.25: 'orange', 0.5: 'green', 0.75: 'red', 1.0: 'purple'}
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        # F = X[-1]
        # Pi = -np.diff(F)/np.diff(D)
        # D = 0.5*(D[1:] + D[:-1])
        # plt.plot(X[0], X[3], 'o')
        plt.scatter(D[1:-1:2], L[1:-1:2],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        if chi == 0.0:
            plt.plot(2*tb.d, [tb.overlap_zone(0.435, d) for d in tb.d], '--', color='black')
    plt.loglog()
    # plt.ylim(1e-3, 5)
    plt.legend()
    plt.savefig('L_on_D.jpg')
    plt.close()
    
    
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        G = X[2]
        # F = X[-1]
        # Pi = -np.diff(F)/np.diff(D)
        # D = 0.5*(D[1:] + D[:-1])
        # plt.plot(X[0], X[3], 'o')
        plt.scatter(D[1:-1:2], G[1:-1:2],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot(2*tb.d, [tb.overlap_integral(0.2, d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    # plt.ylim(1e-3, 5)
    plt.legend()
    plt.savefig('Gamma_on_D.jpg')
    plt.close()
    
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        G = X[2]
        F = X[-1]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        G = 0.5*(G[1:] + G[:-1])
        # plt.plot(X[0], X[3], 'o')
        # plt.scatter(D[0:-1:5], G[0:-1:5]/Pi[0:-1:5],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        # plt.plot(2*tb.d, [tb.overlap_integral(0.2, d)/tb.pressure(d) for d in tb.d], '--', color=colors[chi])
        plt.scatter(Pi[0:-1:5], G[0:-1:5]/Pi[0:-1:5],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot([tb.pressure(d) for d in tb.d], [tb.overlap_integral(0.2, d)/tb.pressure(d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    # plt.ylim(0, 2)
    plt.legend()
    plt.savefig('mu_on_D.jpg')
    plt.close()
    