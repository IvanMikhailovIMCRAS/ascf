from ascf import FlatBrush, TwoBrush
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')

if __name__ == "__main__":
    N = 500
    sigma = 0.1
    colors = {0.0: 'blue', 0.25: 'orange', 0.5: 'green', 0.75: 'red', 1.0: 'purple'}
    
    # L on D
    # plt.figure(figsize=(8.2, 6.0))
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        plt.scatter(D[1:-1:2], L[1:-1:2],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        if chi == 0.0:
            plt.plot(2*tb.d, [tb.overlap_zone(0.435, d) for d in tb.d], '--', color='black')
    plt.loglog()
    plt.yticks([4.0,5.0,6.0],minor=True)
    # plt.ylim(3, 6)
    plt.xlabel('$D$', fontsize=14)
    plt.ylabel('$L(D)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('L_on_D.jpg')
    plt.close()
    
    
    # phi_m on D
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        phi_m = X[3]
        plt.scatter(D[1:-1:5], phi_m[1:-1:5],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot(2*tb.d, [tb.phi_midplace(d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    # plt.ylim(1, 10)
    plt.xlabel('$D$', fontsize=14)
    plt.ylabel(r'$\varphi_m(D)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('phi_m_on_D.jpg')
    plt.close()
    
    # Gamma on D
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        G = X[2]
        plt.scatter(D[1:-1:2], G[1:-1:2],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot(2*tb.d, [tb.overlap_integral(0.2, d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    plt.ylim(1e-3, 5)
    plt.xlabel('$D$', fontsize=14)
    plt.ylabel('$\Gamma(D)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('Gamma_on_D.jpg')
    plt.close()
    
    # Pi on D
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        G = X[-1]
        F = X[-2]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        G = 0.5*(G[1:] + G[:-1])
        plt.scatter(D[0:-1:5], Pi[0:-1:5],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot(2*tb.d, [tb.pressure(d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    plt.ylim(1e-3, 3)
    # plt.xlim(1e-3, 3)
    plt.xlabel('$D$', fontsize=14)
    plt.ylabel('$\Pi(D)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('Pi_on_D.jpg')
    plt.close()
    
    
    # f on Pi
    k = (np.pi - 2) / 4
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        G = X[-1]
        F = X[-2]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        G = 0.5*(G[1:] + G[:-1])
        plt.scatter(Pi[0:-1:5], G[0:-1:5],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot([tb.pressure(d) for d in tb.d], [k*tb.overlap_integral(0.2, d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    plt.ylim(1e-3, 3)
    plt.xlim(1e-3, 3)
    plt.xlabel('$\Pi$', fontsize=14)
    plt.ylabel('$f(\Pi)/(V \zeta)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('f_on_Pi.jpg')
    plt.close()
    
    # mu on Pi
    k = (np.pi - 2) / 4
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N, sigma, chi)
        X = np.loadtxt(f'chi{chi}.txt')
        D = X[0]
        L = X[1]
        G = X[-1]
        F = X[-2]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        G = 0.5*(G[1:] + G[:-1])
        plt.scatter(Pi[0:-1:5], G[0:-1:5]/Pi[0:-1:5],s=12,marker='o', color=colors[chi], label=f"$\chi={chi}$")
        plt.plot([tb.pressure(d) for d in tb.d], [k*tb.overlap_integral(0.2, d)/tb.pressure(d) for d in tb.d], '--', color=colors[chi])
    plt.loglog()
    plt.ylim(0.05, 10)
    plt.xlim(1e-3, 3)
    plt.xlabel('$\Pi$', fontsize=14)
    plt.ylabel('$\mu/(V \zeta)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('mu_on_Pi.jpg')
    plt.close()
    