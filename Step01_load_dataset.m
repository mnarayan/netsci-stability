startup;
tmp_dataset = load('Data/RestingMatrices.mat')
Sigma = tmp_dataset.Sighat; % Specify your own pxp empirical covariance matrix or coherence matrix
usePlotly = 1;
graphtype = 'Resting_p18'
