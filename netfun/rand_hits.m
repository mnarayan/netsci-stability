% rand_hits.m
% function [hubs authorities] = rand_hits(A,varargin)
%  This function does to HITS what pagerank does to katz/eigenvector centrality
%
% See "Stable Algorithms for Link Analysis" by Ng, Zheng and Jordan
% Manjari Narayan
% Copyright 2016
function [hubs authorities] = rand_hits(A,varargin)

  if(nargin>=2)
		opts = varargin{1};
		if(isfield(opts,'visible'))
			visible = opts.visible;
		else
			visible = 0;
		end

		if(isfield(opts,'beta'))
			beta = opts.beta;
		else
			beta = .01;
		end
	else
		visible = 0;
		bipartite = 1;
		beta = .01;
	end

assert(sum(sum(abs(A-A')))>1,'Adjacency matrix is symmetric and undirected. Use eigenvector/katz centrality instead');

p = size(A,1);
h_k = ones(1,p);
a_k = ones(1,p);
u = ones(1,p);

D_row = sum(A,2); D_row(D_row==0) = 1;
D_col = sum(A,1);  D_col(D_col ==0) = 1;

A_row = inv(diag(D_row))*(A);
if(visible)
  disp('Check row stochastic')
  sum(A_row,2)'
end
A_col = (A)*inv(diag(D_col));
if(visible)
  disp('Check column stochastic')
  sum(A_col,1)
end

hasConverged = 0;
iter = 0; max_iter = 1000;

while (hasConverged==0)

  update_h = beta*u + (1-beta)*(A_col*a_k')';
  update_a = beta*u + (1-beta)*(A_row'*h_k')';

 if(sum(abs(h_k-update_h).^2) < 1e-10 & sum(abs(a_k-update_a).^2) < 1e-10)
   if(visible)
     disp(['Converged in ' num2str(iter) ' iterations.'])
   end
    hasConverged=1;
  elseif(iter>max_iter)
    hasConverged = 1;
 end

  h_k = reshape(update_h, [1 p]);
  a_k = reshape(update_a, [1 p]);
  iter = iter+1;
end


hubs = h_k'/sum(h_k);
authorities = a_k'/sum(a_k);
