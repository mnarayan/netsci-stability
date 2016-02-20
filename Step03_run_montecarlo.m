%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Correlation Network %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 20*p;
[V D] = eig(Sigma);
sqrtSig = sqrtm(Sigma);

X = randn(t,p)*sqrt(D)*V';
disp('Verifying simulated Frobenius norm error at t=1000')
frob_err = sum(sum((cov(X)-Sigma).^2,1),2)
assert(frob_err/p<=1,'Error too large')

% X = randn(t,p)*sqrtSig;
% disp('Verifying simulated Frobenius norm error at t=1000')
% frob_err = sum(sum((cov(X)-Sigma).^2,1),2)
% assert(frob_err/p<=1,'Error too large')

X = randn(2*t,p)*sqrt(D)*V';
disp('Verifying simulated Frobenius norm error at t=2000')
frob_err = sum(sum((cov(X)-Sigma).^2,1),2)
assert(frob_err/p<=1,'Error too large')



netstat_global = zeros(n_trials,length(obs),n_thresh);
netstat_nodes = zeros(n_trials,length(obs),p,n_thresh);


if(isWeighted)
	for trial=1:n_trials
		trial
		for t=1:length(obs)
			X = randn(obs(t),p)*sqrt(D)*V';
			Sighat = cov(X); Sighat(find(eye(p))) = 0;
			for tau=1:length(taus)
				%soft(g, τ ) := sign(g) · (|g| − τ )+
				softSig = (abs(Sighat)-taus(tau));
				threshSig = softSig.*(softSig>=0) + 0.*(softSig<0);
				softthreshSig = sign(Sighat).*threshSig;
				switch func2str(bct_funs{bct_num})
				case 'efficiency_bin'
					tmp_stats = feval(bct_funs{bct_num},softthreshSig,1);
				case 'efficiency_wei'
					tmp_stats = feval(bct_funs{bct_num},softthreshSig,1);
				otherwise
					tmp_stats = feval(bct_funs{bct_num},softthreshSig);
				end
				netstat_nodes(trial,t,:,tau) = tmp_stats;
			end
		end
	end
else
	for trial=1:n_trials
		trial
		for t=1:length(obs)
			X = randn(obs(t),p)*sqrt(D)*V';
			Sighat = cov(X); Sighat(find(eye(p))) = 0;
			for tau=1:length(taus)
				tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));
				netstat_nodes(trial,t,:,tau) = tmp_stats;
			end
		end
	end
end

% Plot all nodes
% for time=1:length(obs)
% 	subplot(length(obs),1,time);
% 	title(sprintf('t = %3.1f',obs(time)),'fontsize',16);
% 	xlabel('Varying threshold');
% 	ylabel('Network Metric')
% 	tmp_x = reshape(taus,[n_thresh 1]);
% 	tmp_y = squeeze(mean(netstat_nodes(:,time,:,:),1)); tmp_y = reshape(tmp_y,[n_thresh p]);
% 	plot(repmat(tmp_x,[1 p]),tmp_y);
% 	ylim([0 25]);
% 	xlim([.08 .22]);
% end
