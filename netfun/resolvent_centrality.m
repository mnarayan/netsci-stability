% resolvent_centrality.m
% function centrality = resolvent_centrality(A,varargin)
% INPUTS
% A is an adjacency matrix for a weighted or unweighted graph
% optional:
% opts.beta is the dampening factor (regularization factor)
%
% M = (I-beta*A)^-1 * d
%   where A is a stochastic adjacency matrix
%   			d is the degree vector
% 				beta is the dampening factor. Can be positive or negative
% Thus dominant eigenvector of M is the bonanich power centrality
% 
% When abs(beta) is 1/max_eigenvalue(A) and D is the all ones vector
% the dominant eigenvector of M is the eigenvector centrality
% 
% When 0 <= beta <= 1/max_eigenvalue(A) and D is the all ones vector, 
% the dominant eigenvector of M is the katz centrality
% 
% Manjari Narayan
% Copyright 2017

function centrality = resolvent_centrality(A,varargin)

  if(nargin>=2)
		opts = varargin{1};
		if(isfield(opts,'visible'))
			visible = opts.visible;
		else
			visible = 0;
		end

		if(isfield(opts,'damping'))
			damping_factor = opts.damping;
		else
			damping_factor = .99;
		end
		
		if(isfield(opts,'resolvent_type'))
			resolvent_type = opts.resolvent_type;
		else
			resolvent_type = 'eigenvector';
		end
	else
		visible = 0;
		damping_factor = .99;
		resolvent_type = 'eigenvector'; 
	end

	%assert(sum(sum(abs(A-A')))< 1e-3,'Adjacency matrix is asymmetric and directed. Use pagerank or rand_hits instead');

	if(visible)
		disp('Forcing A to be unsigned')
	end
	p = size(A,1);
	A = abs(A); 
	D_col = sum(abs(A),2)';  D_col(D_col ==0) = 1;
	Deig = eigs(A); 
	max_eigenvalue = max(Deig); 


	% if(beta==0)
	% 	beta = .99/max_eigenvalue;
	% end

	% update_centrality = zeros(1,p);
	% curr_centrality = zeros(1,p);
	%
	% hasConverged = 0;
	% iter = 0; max_iter = 1000;
	%
	% while (hasConverged==0)
	%   update_centrality = D_col + (beta*A*curr_centrality')';
	%  if(sum(abs(curr_centrality-update_centrality).^2) < 1e-10)
	%    if(visible)
	%      disp(['Converged in ' num2str(iter) ' iterations.'])
	%    end
	%     hasConverged=1;
	%  elseif(iter>max_iter)
	%   hasConverged = 1;
	%  end
	%   curr_centrality = reshape(update_centrality, [1 p]);
	%   iter = iter+1;
	% end
	% centrality = update_centrality'/sum(update_centrality);

	% Calculate matrix resolvent using linsolve(A,b,opts)
	% AX = b
	% 
	linopts = {};
	linopts.SYM = true;	
	switch resolvent_type
	case 'eigenvector'
		% Standard eigenvector
		if(visible)
			disp('Standard eigenvector centrality')
		end
		damping_factor = .99;
		beta = abs(damping_factor)/max_eigenvalue;
		M = eye(p)-beta*A;
		b = ones(p,1);
	case 'katz'
		% Katz type
		%damping_factor = .8;
		disp('Katz centrality')
		beta = abs(damping_factor)/max_eigenvalue;
		M = eye(p)-beta*A;
		b = ones(p,1);
	case 'power'
		% Bonanich Power
		disp('Bonanich power centrality')
		%damping_factor = .8;
		beta = damping_factor/max_eigenvalue;
		M = eye(p)-beta*A;
		b = A*ones(p,1);
	case 'walk_probability'
		%	Random Walk
		disp('Random Walk (pagerank) centrality')
		invD = inv(diag(D_col));
		M = eye(p)*(1-abs(damping_factor)) + abs(damping_factor)*(eye(p)-A*invD);
		b = ones(p,1);
	otherwise
	% Standard eigenvector
		disp('Standard eigenvector centrality')
		damping_factor = .99;
		beta = abs(damping_factor)/max_eigenvalue;
		M = eye(p)-beta*A;
		b = ones(p,1);
	end
	
	centrality = linsolve(M,b); 
	normalization = 1.0; 
	% normalization = sum(centrality); 
	centrality = centrality'/normalization; 
	if(visible)
		disp('Top central nodes')
		find(centrality>=.99*max(centrality))
	end

end
