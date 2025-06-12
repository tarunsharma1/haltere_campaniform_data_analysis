import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import pickle

def estimate_amplitude(sig, dt, frequency, num_cycle, t):
    sin_ref = np.sin(2.0*np.pi*frequency*t)
    cos_ref = np.cos(2.0*np.pi*frequency*t)
    int_sig_sqr = sp.integrate.trapezoid((sig - sig.mean())**2, dx=dt)
    return np.sqrt(2*frequency*int_sig_sqr/num_cycle)


def estimate_phase(sig, dt, frequency, num_cycle, t):
    sin_ref = np.sin(2.0*np.pi*frequency*t)
    cos_ref = np.cos(2.0*np.pi*frequency*t)
    sig_sub_mean = sig - sig.mean()
    int_sig_cos = sp.integrate.trapezoid(sig_sub_mean*cos_ref, dx=dt)
    int_sig_sin = sp.integrate.trapezoid(sig_sub_mean*sin_ref, dx=dt)
    return np.arctan2(int_sig_cos, int_sig_sin)
    

# ------------------------------------------------------------------------------
if __name__ == '__main__':

    num_pts = 1000
    num_cycle = 3
    frequency = 1.0
    amplitude = 2.0
    offset = 1.2
    phase_deg = 65.0
    phase_rad = np.deg2rad(phase_deg)
    
    period = 1/frequency
    duration = num_cycle*period
    
    t = np.linspace(0.0, duration, num_pts)
    
    #sig = amplitude*np.cos(2.0*np.pi*frequency*t + phase_rad) + offset
    
    t,sig = pickle.load(open('/home/tarun/Desktop/curve1.p', 'rb'))
    sig = sig - np.mean(sig)
    dt = t[1] - t[0]
    frequency = 1/(18.01)
    num_cycle = 0.5


    amplitude_est = estimate_amplitude(sig, dt, frequency, num_cycle, t)
    phase_rad_est = estimate_phase(sig, dt, frequency, num_cycle, t)
    phase_deg_est = np.rad2deg(phase_rad_est)
    

    sig_reconstructed = amplitude_est*np.sin(2.0*np.pi*frequency*t + phase_rad_est)
    plt.plot(t, sig, 'r')
    plt.plot(t, sig_reconstructed, 'b')
    plt.show()



    print()
    print(f'amplitude:     {amplitude}')
    print(f'amplitude_est: {amplitude_est}')
    print()
    print(f'phase_rad:     {phase_rad}')
    print(f'phase_rad_est: {phase_rad_est}')
    print()
    print(f'phase_deg:     {phase_deg}')
    print(f'phase_deg_est: {(phase_deg_est - 90)%360}')
    print()














