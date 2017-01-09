% centrality2centralization
% function centralization = centrality2centralization(centrality, centrality_type)
% Centralization is a global network level metric. In contrast to global network metrics that capture the mean node level metric by taking the average centrality value across all nodes, centralization captures variation in centrality across all nodes in the network. For a given centrality metric, centralization measures whether all nodes share similar values of centrality (centralization towards 0) or whether certain nodes have large centrality while others have lower centrality (centralization towards 1). 
% 
% INPUTS: 
% centrality 	- a nodes x 1 length vector of centrality values per node in the graph
% centrality_type 	- Currently supports 'degree', 'betweenness', 'closeness', 'eigenvector'. General support coming. Is not defined for valued networks except for information centrality or current flow closeness centrality. 
% normalize 	- If 1 normalize by theoretical maximum centralization. Otherwise return unnormalized

% To Be Done: Centralization for any centrality can be implemented using a wrapper. 
% 1. Call centrality2centralization for given network and centrality function, 
% 2. Then call centrality2centralization for theoretically most central network i.e. star topology of size p
% 3. Divide centralization from Step 1 / centralization from Step 2

function [centralization diff_centrality] = centrality2centralization(centrality, centrality_type, normalize, varargin)

	
if(length(size(centrality))>2)
	disp('Input should be a vector of length (no. of nodes x 1)');
end
	
if(nargin==4)
	cfun = varargin{1}; 
else
	cfun = [];
end
	
p = length(centrality)	;
max_centralization = 0;
	
most_central = max(centrality); 
diff_centrality = abs(bsxfun(@minus, most_central, centrality));


% TBD: Create star graph and compute maximum possible centrality scores for the star graph. 
if(normalize)
	max_G = zeros(p,p); max_G(1,2:end) = 1; max_G(2:end,1) = 1;

	switch centrality_type
	
	case 'degree'
		tmp_metric = sum(abs(max_G)); 
		max_centralization = sum(bsxfun(@minus,max(tmp_metric),tmp_metric));
	
	case 'betweenness'
		% maximum occurs for a star topology
	 	max_centralization = (p-1)*nchoosek(p-1,2);
	case 'closeness'
		% maximum occurs for a star topology
		max_centralization = (p-2)*(p-1)/(2*p-3);	
	case 'eigenvector'
		max_G(find(eye(p))) = sum(max_G);
		[tmpV tmpD] = eig(max_G);
		[max_val max_idx] = max(diag(tmpD));
		tmp_metric = tmpV(:,max_idx);
		max_centralization = sum(bsxfun(@minus,max(tmp_metric),tmp_metric));
	case 'hub'
		max_G = triu(max_G,1);
		[hubs authorities eig_centrality] =  hits(max_G);
		tmp_metric = hubs;
		max_centralization = sum(bsxfun(@minus,max(tmp_metric),tmp_metric));
		
	case 'authority'
		max_G = tril(max_G,-1);
		[hubs authorities eig_centrality] =  hits(max_G);
		tmp_metric = authorities;
		max_centralization = sum(bsxfun(@minus,max(tmp_metric),tmp_metric));
		
	case 'pagerank'
		max_G = triu(max_G,1);
		beta = .95;
		tmp_metric = pagerank_centrality(max_G, beta);
		max_centralization = sum(bsxfun(@minus,max(tmp_metric),tmp_metric));		
		
	case 'custom'
		disp('Not implemented')
			
	end
else
	max_centralization = (p-1)*max(diff_centrality); %norm(diff_centrality); %most_central;
	if(max_centralization==0)
		max_centralization = 1.0;
	end
end

centralization = sum(diff_centrality)/max_centralization;
	
end