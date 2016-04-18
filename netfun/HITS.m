function [hubs authorities eig_centrality] =  HITS(A,varargin); 
	% HITS Algorithm 
	% (Kleinberg)
	% Given a directed and thus asymmetric connectivity matrix A, HITS(A) returns centrality scores for hubs and authorities
	% These are directed analogues of eigenvector centrality
	% 
	% 
	% Test matrix
	% A = zeros(4,4);
	% A(1,2) = 1; A(1,3) = 1;
	% A(3,4) = 1; A(3,2) = 1;
	% A(4,1) = 1;
	%
	%
	% Test matrix
	% A = zeros(3,3); 
	% A(1,3) = 1; A(2,3) = 1; 
	% hubs = (2,2,0); authorities = (0,0,2);

	if(nargin>=2)
		opts = varargin{1};
		if(isfield(opts,'visible'))
			visible = opts.visible;
		else
			visible = 0;
		end
		
		if(isfield(opts,'bipartite'))
			bipartite = opts.bipartite;
		else
			bipartite = 1;
		end
		
		if(isfield(opts,'beta'))
			beta = opts.beta;
		else
			beta = 1;
		end
	else
		visible = 0;
		bipartite = 1;
		beta = 1;
	end
	
	assert(sum(sum(abs(A-A')))>1,'Adjacency matrix is symmetric and undirected');
	
	L_in = beta.*A*A' + (1-beta)/size(A,1).*ones(size(A,1),size(A,1));
	L_out = beta.*A'*A + (1-beta)/size(A,1).*ones(size(A,1),size(A,1));
	
	[V_in D_in] = eig(L_in); 	
	[sorted_D_in idx_in]  = sort(diag(D_in),'descend');
	if(visible)
		disp('Largest two eigenvalues of AA^T ')
		sorted_D_in(1:2)'	 
	end
	hubs = V_in(:,idx_in(1)); % ordered from highest to lowest
	hubs = hubs/sum(hubs);
	assert(D_in(idx_in(1),idx_in(1))~=D_in(idx_in(2),idx_in(2)),'Principal eigenvalue is not unique. No stable hub/auth score');
	
	[V_out D_out] = eig(L_out); 	
	[sorted_D_out idx_out]  = sort(diag(D_out),'descend'); 
	if(visible)
		disp('Largest two eigenvalues of A^TA ')
		sorted_D_out(1:2)'
	end
	authorities = V_out(:,idx_out(1)); % ordered from highest to lowest
	authorities = authorities/sum(authorities);
	
	% Make symmetric, 
	symA = (A+A'); 
	diagA = inv(diag(sum(abs(symA))));
	symA = diagA.^(.5)*symA*diagA.^(.5);
	[V D] = eig(symA);
	[sorted_D idx_und]  = sort(diag(D),'descend');
	if(visible)
		disp('Largest two eigenvalues of symmetricized A ')
		sorted_D(1:2)'	 		
	end
	eig_centrality = V(:,idx_in(1)); % ordered from highest to lowest
	eig_centrality = eig_centrality /sum(eig_centrality );

	