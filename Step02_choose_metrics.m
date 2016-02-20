%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick Metrics %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
bct_funs = {@betweenness_bin, ... %1
									@clustering_coef_bu, ... 	%2
									@efficiency_bin, ...
									@eigenvector_centrality_und, ...
									@betweenness_wei, ...
									@clustering_coef_wu, ...
									@efficiency_wei, ...
									 };
bct_num = 7;
isWeighted = 0; % Eventually choose this and network metrics with error checking
disp(['Metric Chosen is ' func2str(bct_funs{bct_num}) ', isWeighted: ' num2str(isWeighted)])


%%%%%%%%%%%%%%%%%%
%%%Choose Figure Parameters%%%
%%%%%%%%%%%%%%%%%%%
p = size(Sigma,1);
node = floor(p/2); % Specity node/region of your choice
t_idx = 5; % Choose from 1 to 5 where 1 corresponds to few observations at t=2p, and 5 corresponds to t=10p



%%%%%%%%%%%%%%%%%%
%%%Choose Simulation Parameters%%%
%%%%%%%%%%%%%%%%%%%
obs = round(linspace(2*p,10*p,5)); % sample sizes: 5 values spaced between 2p and 10p
%%%%%%%%%%%%%%%%%%%
n_trials = 50; % Monte carlo trials
%%%%%%%%%%%%%%%%%%%
n_thresh = 25; % Granularity of shrinkage or thresholding parameter
tau_start = .1;
tau_stop = .2;
taus = linspace(tau_start,tau_stop,n_thresh);
%%%%%%%%%%%%%%%%%%%
