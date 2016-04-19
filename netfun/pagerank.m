% pagerank.m
% function centrality = pagerank(A,varargin)
% INPUTS
% A is an adjacency matrix for a directed graph
% optional:
% opts.beta is the dampening factor (regularization factor)
%
% M = (beta*U + (1-beta)*A_row)'
%   where A_row is a row stochastic adjacency matrix
%   where U is the all ones matrix normalized by 1/(no. of nodes)
% Thus eigenvector of M is the pagerank centrality
%
% See "Stable Algorithms for Link Analysis" by Ng, Zheng and Jordan
% Manjari Narayan
% Copyright 2016
function centrality = pagerank(A,varargin)

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
pagerank = ones(1,p);
u = ones(1,p);

D_col = sum(A,1);  D_col(D_col ==0) = 1;

A_col = (A)*inv(diag(D_col));
if(visible)
  disp('Check column stochastic')
  sum(A_col,1)
end

hasConverged = 0;
iter = 0; max_iter = 1000;

while (hasConverged==0)

  update_pg = beta*u + (1-beta)*(A_col*pagerank')';

 if(sum(abs(pagerank-update_pg).^2) < 1e-10)
   if(visible)
     disp(['Converged in ' num2str(iter) ' iterations.'])
   end
    hasConverged=1;
elseif(iter>max_iter)
  hasConverged = 1;
 end

  pagerank = reshape(update_pg, [1 p]);
  iter = iter+1;
end


centrality = pagerank'/sum(pagerank);


end
