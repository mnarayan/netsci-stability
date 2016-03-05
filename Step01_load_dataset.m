startup;

datapath = 'Data'
datafile = 'RestingMatrices.mat'
fieldname = 'Sighat' % Variable containing an empirical correlation, covariance matrix;
useParCor = 0; % Using standard correlation,covariance matrix

try
	help plotly;
	usePlotly = 1; 	% Use plotly toolbox if available
catch
	usePlotly = 0;	% Plotly likely unavailable.
end

visible = 0;	% Plot figures or turn off plot checks.

if(~useParCor) 
	graphtype = 'Resting_p18'; % Name dataset/graph origin. 
else
	graphtype = 'Resting_p18_parcor'; % Name dataset/graph origin. 
end

tmp_dataset = load([datapath '/' datafile]);
Sigma = getfield(tmp_dataset,fieldname); % Specify your own pxp empirical covariance matrix or coherence matrix
