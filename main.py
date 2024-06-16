from ascf import FlatBrush, TwoBrush
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    # for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
    #     fb = FlatBrush(N=100, sigma=0.1, chi=chi)
    #     plt.plot(fb.z, fb.phi, label=f"$\chi={chi}$")
    # plt.legend()
    # plt.show()
    
    for chi in [0.0, 0.25, 0.5, 0.75, 1.0]:
        tb = TwoBrush(N=100, sigma=0.1, chi=chi)
        # plt.plot(2*tb.d, [tb.phi_midplace(d) for d in tb.d], label=f"$\chi={chi}$")
        plt.plot(2*tb.d, [tb.pressure(d) for d in tb.d], label=f"$\chi={chi}$")
    plt.loglog()
    plt.legend()
    plt.show()
    