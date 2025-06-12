Code for analyzing head and wing equilibrium responses in Drosophila, across various genotypes that perturb increasing number of campaniform sensilla using either Kir2.1 or reaper,hid

Data was collected by subjecting flies to a angular velocity stimulus about their yaw axis. A rotating beam, connected to a servo motor and controlled using the GRBL library in python, consisting of a tether holder, IR light source, camera, and IMU was used for these experiments. Flies were subjected to 5 sinusoidal cycles reaching a maximum angular velocity of 350 deg/sec, followed by a 15 second rest period and a second set of 5 sinusoidal cycles reaching a maximum angular velocity of 500 deg/sec. Experiments were performed on 9 screened GAL4 lines targeting various subsets of campaniform sensilla on the halteres. Lines were crossed with Kir2.1 and reaper-hid, to test the effects of silencing or ablation of campaniforms on the wing and head equilibrium responses. Data from all sensors was stored in bag files using ROS. Left and right wing beat amplitudes were measured using Kinefly and head yaw movements were tracked using DLC by tracking points on the edges of the head. Three trials were performedper fly and data averaged. Data from n=18-20 flies per genotype were collected. 

bag_to_hdf5.py: converts bag files into compressed hdf5 files for ease of data analysis.

rkf_analysis_tarun.py: is run for individual flies. For each trial per fly, the code segregates the trajectory, wing, and head data into the separate set of 5 cycles each, applies a butterworth filter, and resamples wing, head, and trajectory data to a common frequency. The three trials per fly are then resampled to a common time vector, averaged and stored in pickle files.

rkf_analysis_tarun_plotting_all_flies.py: this script reads the saved pickle files for all flies of a genotype (one file fly averaged over 3 trials), resamples them to a common time stamp, averages and then plots everything on one plot. For each fly and each set of 5 sinusoidal cycles, data is collapsed over the 5 cycles (by finding the points of zero crossings of the trajectory), and the wing and head responses averaged over the 5 cycles are quantified using a fourier transform. The amplitude and phase of the fourier transform fits are saved on a per fly basis for comparisons between populations (genotypes).

fft_check.py: Test script for visualizing the fits using fourier transform to a noisy signal. 
will_amplitude_and_phase_calc.py: Alternate way for calculating fits to a noisy signal.

test_for_normality.py and mann_whitney_u_test.py: Stastical tests comparing responses from populations (genotypes).

inside_outside_wing_comparisons.py: Testing the effects of silencing the campaniforms on the inside wing (increasing in amplitude) and outside wing (decreasing in amplitude).

least_squares_solution_campaniforms.py: Least squares regression to determine the line of best fit to wing or head responses plotted against number of campaniforms silenced.

linear_regression_across_fields.py: Regression analysis to determine the effect of any individual field on linear trend observed.

scatter_plot.py: Code to visualize amplitudes of wing or head responses across different genotypes using a scatter plot.

wing_vs_head_plots.py and wing_vs_head_subplots.py: Code to visualize wing responses vs head responses.

lissajous.py: Code to make lissajous figure of wing vs head phase offset.

 

