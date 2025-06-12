import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.signal as sps
from scipy.optimize import leastsq

def fourier_transform_to_get_cutoff(signal, t):

	### one alternative suggested by email to reduce noise was to take the second half of the signal,invert it as the first half and use this as the cleaned signal
	# N = signal.shape[0]
	# signal_second_half = signal[int(N/2)::]
	# ## move the signal down to 0 to fix asymmetry
	# signal_second_half = signal_second_half - signal_second_half[0]
	# signal_new = np.zeros_like(signal)
	# signal_new[0:(int(N%2) + int(N/2))] = -1.0*signal_second_half
	# signal_new[int(N/2)::] = signal_second_half

	# signal = signal_new



	# helper function for me to visualize and see 
	# fourier transform to freq domain and look at signal to see what cutoff value to use for low pass filter

	N = signal.shape[0]
	dt = t[2] - t[1]

	#print ('sampling freq:', 1.0/dt)

	fft = 1.0/N * np.fft.fft(signal)
	fft = fft[:N//2]
	fftfreq = np.fft.fftfreq(N, dt)[:N//2] 
	
	#plt.plot(fftfreq, np.abs(fft))
	#plt.show()
	
	# the first highest fft value is going to be a line ..look at fun with fft if you set the signal as 5 + sinx biggest spike is a line then is the one for the sine curve
	# so we can safely ignore the first fft value.. and we know its going to be at 0 freq i.e at fft[0]

	idx = 1 + np.argmax(np.abs(fft[1:]))
	#print ('index:', idx)

	amp = np.abs(fft[idx]) * 2
	#print ('amp:', amp)

	phase =  np.angle(fft[idx]) * 180/np.pi
	#print ('phase in degrees:',phase)

	freq = np.abs(fftfreq[idx])
	#print ('frequency:',freq)

	## add offset np.mean(signal) just for visualizing the fits better
	#phase = 90
	signal_reconstucted =  amp*np.cos(2*np.pi*freq*t + phase*np.pi/180) + np.mean(signal)


	### lets try doing an optimization here using scipy to increase fit quality
	optimize_func = lambda x: x[0]*np.cos(2*np.pi*x[1]*t+x[2]*np.pi/180) + np.mean(signal) - signal
	est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [amp, freq, phase, np.mean(signal)])[0]
	#print ('##### scipy optimized #####')
	#print ('amp:', est_amp)
	#print ('phase', est_phase)
	#print ('frequency:', est_freq)
	optimized_signal_reconstructed =  est_amp*np.cos(2*np.pi*est_freq*t + est_phase*np.pi/180) + np.mean(signal)

	plt.clf()
	#plt.ylim(-8,8)
	return optimized_signal_reconstructed, est_amp, est_phase, signal
	#return signal_reconstucted, amp, phase, signal


#t = np.arange(0,2*np.pi, 1/72)
#curve = -1.5*np.sin(t) + np.random.normal(0,0.2,t.shape[0])

#curve_half = curve[int(curve.shape[0]/2)::]
#t = t[int(curve.shape[0]/2)::]

t,curve = pickle.load(open('/home/tarun/Desktop/curve1.p', 'rb'))


signal_reconstucted, amp, phase, signal = fourier_transform_to_get_cutoff(curve, t)
print (f'amp: {amp} and phase is {phase}')
plt.plot(t, signal, 'r')
plt.plot(t, signal_reconstucted, 'b')
plt.show()