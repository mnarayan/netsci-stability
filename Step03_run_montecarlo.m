%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Correlation Network %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 20*p;
[V D] = eig(Sigma);

Z = randn(t/2,p);
X = Z*V*sqrt(D)*V';
Y = Z*sqrtm(Sigma);
disp(['Verifying simulated Frobenius norm error at t=' num2str(t/2)])
frob_err = sum(sum((cov(X)-Sigma).^2,1),2)
assert(frob_err/p^2<=.1,'Error too large')
assert(sum(sum((cov(X)-cov(Y)).^2,1),2)<=eps,'Simulation Incorrect')


X = randn(t,p)*V*sqrt(D)*V';
disp(['Verifying simulated Frobenius norm error at t=' num2str(t)])
frob_err = sum(sum((cov(X)-Sigma).^2,1),2)
assert(frob_err/p^2<=1e-2,'Error too large')



netstat_global = zeros(n_trials,length(obs),n_thresh);
netstat_nodes = zeros(n_trials,length(obs),p,n_thresh);


if(isWeighted)
	for trial=1:n_trials
		trial
		for t=1:length(obs)
			X = randn(obs(t),p)*V*sqrt(D)*V';
			Sighat = cov(X); Sighat(find(eye(p))) = 0;
			for tau=1:length(taus)
				%soft(g, τ ) := sign(g) · (|g| − τ )+
				softSig = triu((abs(Sighat)-taus(tau)),1);
				threshSig = softSig.*(softSig>=0) + 0.*(softSig<0);
				softthreshSig = sign(Sighat).*threshSig;
				softthreshSig = softthreshSig + softthreshSig';
				assert(sum(diag(softthreshSig)<=0)~=0,'Negative values on diagonal');
				%[V2 D2] = eig(softthreshSig); softthreshSig = V2*(D2+eye(p)*(.01+abs(min(diag(D2)))))*V2';
				switch func2str(bct_funs{bct_num})
				case 'efficiency_wei'
					tmp_stats = feval(bct_funs{bct_num},softthreshSig,1);
				otherwise
					tmp_stats = feval(bct_funs{bct_num},softthreshSig);
				end
				if(~isreal(tmp_stats))
					warning('Network metrics are complex. Taking absolute values'); 
					tmp_stats = abs(tmp_stats);
				end
				netstat_nodes(trial,t,:,tau) = tmp_stats;
			end
		end
	end
else
	for trial=1:n_trials
		trial
		for t=1:length(obs)
			X = randn(obs(t),p)*V*sqrt(D)*V';
			Sighat = cov(X); Sighat(find(eye(p))) = 0;
			assert(sum(diag(Sighat)<=0)~=0,'Negative values on diagonal');			
			for tau=1:length(taus)
				switch func2str(bct_funs{bct_num})
				case 'efficiency_bin'
					tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)),1);
				otherwise
					tmp_stats = feval(bct_funs{bct_num},1*(abs(Sighat)>taus(tau)));
				end
				if(~isreal(tmp_stats))
					warning('Network metrics are complex. Taking absolute values'); 
					tmp_stats = abs(tmp_stats);
				end
				netstat_nodes(trial,t,:,tau) = tmp_stats;
			end
		end
	end
end

save('tmp_netstat_nodes','netstat_nodes');

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
