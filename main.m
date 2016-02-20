% continuity_covariance_networks.m
%
% Questions relevant to regular correlation/partial correlation matrices.
% Q1: Plot network metrics as a function of hard thresholding parameter & binary connections
% 		-
% Q2: Plot sampling variability of network metric as a function of the number of samples t
%
% Q3: Plot network metrics for fully weighted matrics as a function of soft thresholding parameter

Step01_load_dataset;

Step02_choose_metrics;

Step03_run_montecarlo;

set(0, 'DefaultFigureVisible','off');

Step04_plot_metrics;

Step05_cleanup;
