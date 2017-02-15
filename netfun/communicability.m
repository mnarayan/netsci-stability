classdef communicability
% All communicability related metrics
% Manjari Narayan
% Copyright 2017

	methods(Static)
	
		function [output varargout] = communicability_matrix(A,varargin)
			% communicability_matrix.m
			% function output = communicability(A,varargin)
			% 
			% communicability(u,v) = [exp^A]_(u,v)
			% 
			% INPUTS
			% A is an adjacency matrix for a weighted or unweighted graph
			% optional:
			%			% 
			% OUTPUTS:
			% output is a matrix of same dimensions as A
			output = expm(A); 
		
		end
		
		function [output varargout] = subset(A, G, source, target, varargin)
			% communicability.subset.m
			% function output = communicability.subset(A,G, source, target, varargin)
			% 
			% communicability(u,v) = [exp^A]_(u,v)
			% 
			% INPUTS
			% A is an adjacency matrix for a weighted or unweighted graph
			% optional:
			%			 
			% OUTPUTS:
			% output is a matrix of same dimensions as A
			
			if(~isempty(G))
				G = communicability.communicability_matrix(A); 
			end
			
			output = sum(sum(abs(G(source,target))));
		
		end
		
	
		function [output varargout] = distance(A,G,varargin)
			% communicability.distance (Communicability Distance)
			% function output = communicability.distance(A,varargin)
			% 
			% communicability_distance(u,v) = 
			% 
			% INPUTS
			% A is an adjacency matrix for a weighted or unweighted graph
			% optional:
			%		
			% OUTPUTS:
			
			if(isempty(G))
				G = communicability.communicability_matrix(A);
			end				 
			p = size(G,1); 			
			D = diag(G);
			if(sum(D<=0))
				D = 1+D;
			end
			
			normalized=true;
			output = zeros(size(G)); 
			for ii=1:p
				for jj=1:ii
					output(ii,jj) = D(ii) + D(jj) - 2*G(ii,jj); 
					if(normalized)
						output(ii,jj) = output(ii,jj)/(sqrt(D(ii))*sqrt(D(jj)));
					end
				end
			end
			
			output = output + output';
			output(find(eye(p))) = 0; 
			
		end
	
		function [output varargout] = degree(A,G,varargin)
			% total communicability
			
			if(isempty(G))
				G = communicability.communicability_matrix(A);
				G
			end	
			
			G(find(eye(size(G,1)))) = 0;			
			output = sum(G); 
			
			% rescale from 0 to 1
			% output = bsxfun(@minus,max(output),output)./(max(output)-min(output)); 
			
			% rescale by sum
			% output = output/sum(output); 
			
			output = output/max(output);			
		end
	
		% function output = subgraph(A,G,varargin)
		% 	% subgraph centrality
		%
		% 	if(isempty(G))
		% 		G = communicability.communicability_matrix(A);
		% 	end
		% 	output = diag(G)';
		% 	output = output/sum(output);
		%
		% end
		
		function output = neighborhood(A,G, varargin)
			% output = neighborhood(A,G,varargin)
			% INPUTS
			% - A 
			% - G
			% - method (optional): 'subgraph' or straightup 'eigenvector' of communication distance
			
			if(isempty(G))
				G = communicability.communicability_matrix(A);
			end
			
			if(nargin==3)
				method = varargin{1}; 
			else
				method = 'subgraph'; 
			end

			switch method
			case 'subgraph'
				centrality = diag(G)'; 
				centrality = centrality/sum(centrality);
			case 'eigenvector'
				%centrality = eigenvector_centrality_und(G);			
				centrality = bonanich_power_centrality(G)'; 	
			end
			
			output = centrality; 
			
		end
		
		function output = flowmatrix(A,K,node,varargin)
			% F(s,t)i=12∑jAi,j|G+i,s − G+i,t + G+j,t − G+j,s|
			
			p = size(A,1); 
			output = zeros(p,p); 
			
			if(isempty(node)|nargin<3)
				error('node argument missing')
			end
			
			%K(find(eye(size(K,1)))) = 0;		
			not_node = setdiff([1:p],node);
			for ss=1:p
				for tt=ss:p
					try
						res = K(node,ss) - K(node,tt) + K(not_node,tt) - K(not_node,ss);
					catch
						disp('Sizes not matching')
						res = K(node,ss) - K(node,tt) + K(not_node,tt)' - K(not_node,ss)';						
					end
					output(ss,tt) = .5*sum(A(node,not_node).*abs(res)');
				end
			end
			
		end
		
		
		function output = betweenness(A,varargin)

			p = size(A,1); 
			output = {};
			
			switch nargin
				
			case 2
				nodelist = varargin{1}; 
				source = 1:p;
				target = 1:p;
			case 3				
				source = varargin{1};
				target = varargin{2};
				nodelist = 1:p;
			case 4
				source = varargin{1};
				target = varargin{2};
				nodelist = varargin{3}; 
			otherwise
				nodelist = 1:p;
				source = 1:p;
				target = 1:p;
			end	
		
			
			centrality = sparse(zeros(length(A),1));
			dyadic = cell([1 length(nodelist)]);
			normalization = cell([1 length(nodelist)]);			
			if(isempty(nodelist))
				disp('Skipping betweenness')
			else
				G = communicability.communicability_matrix(A); 
				for ii=1:length(nodelist)
					A_ii = A;
					if(isstruct(nodelist)|iscell(nodelist))
						n_nodes = p-length(nodelist{ii});
						C = (length(source)-1)*(length(target)-1) - (n_nodes);
						A_ii(nodelist{ii},:) = 0; A_ii(:,nodelist{ii}) = 0;	
						G_ii = communicability.communicability_matrix(A_ii);
						betweenness = (G-G_ii)./G;
						betweenness(find(eye(p))) = 0;						
						betweenness(nodelist{ii},nodelist{ii}) = 0;
						centrality(nodelist{ii}) = .5*sum(sum(betweenness(source,target)))/C + ...
						 						.5*sum(sum(betweenness(target,source)))/C;	
														
					else
						n_nodes = p-1;
						C = (length(source)-1)*(length(target)-1) - (n_nodes);
						A_ii(nodelist{ii},:) = 0; A_ii(:,nodelist(ii)) = 0;	
						G_ii = communicability.communicability_matrix(A_ii); 
						betweenness = (G-G_ii)./G;
						betweenness(find(eye(p))) = 0;
						betweenness(nodelist(ii),nodelist(ii)) = 0;						
						centrality(nodelist(ii)) = .5*sum(sum(betweenness(source,target)))/C + ...
												.5*sum(sum(betweenness(target,source)))/C;	
					end				
					dyadic{ii} = sparse(betweenness);
					normalization{ii} = C;											
				end
			end
			
			assert(sum(centrality<0)==0,'Negative values of Betweennes. Not possible')

			output.centrality = centrality;
			output.nodelist = nodelist;	
			output.dyadic = dyadic;		
			output.normalization = normalization;
						
		end
		
	
		function [output centrality]= invariant(A,source,target,name,normalized,varargin)
			
			if(~exist('name','var'))
				name = 'closeness';
			end			
			if(~exist('normalized','var'))
				normalized = false;
			end
			
			centrality = zeros(length(A),1);
			reduced = communicability.reduced_source_target(A,source,target);	
			G = communicability.communicability_matrix(A);
			
			switch name
				
			case 'estrada'
				% estrada type
				centrality(reduced.node_idx) = log1p(diag(G(reduced.node_idx,reduced.node_idx)));
				invariant = sum(centrality); 
				if(normalized)
					invariant = invariant/length(reduced.A);
				end
			
			case 'total'
				% total_comm type
				G(find(eye(length(G)))) = 0;
				centrality(reduced.node_idx) = sum(G(reduced.node_idx,reduced.node_idx));
				if(isempty(source) && isempty(target))
					invariant = sum(centrality);
					if(normalized)
						invariant = invariant/((length(reduced.A)-2)*(length(reduced.A)-1)); 
						centrality = centrality/(length(reduced.A)-1);
					end
				else
					st_idx = cat(1,source(:),target(:));					
					centrality(reduced.node_idx) = sum(G(reduced.node_idx,st_idx),2);
					invariant = sum(centrality); 
					if(normalized)
						invariant = invariant/(length(reduced.A)*(length(st_idx)));
						centrality = centrality/(length(st_idx));
					end
				end
			
			case 'connectivity'				
				G(find(eye(length(A)))) = 0;
				connectivity = G(source,target);				
				invariant = sum(sum(connectivity)); 
				if(normalized)
					invariant = invariant/(length(source)*length(target));
				end
				centrality = invariant*ones(length(G),1); % centrality is not defined for this option. 
			
			case 'closeness'
				K = communicability.distance(reduced.A,G(reduced.node_idx,reduced.node_idx)); 
				st_idx = cat(1,reduced.source(:),reduced.target(:));
				centrality(reduced.node_idx) = communicability.closeness(reduced.A,K,st_idx,normalized);  				
				invariant = sum(centrality);
			end
			
			output = invariant;			
			
		end
		
		
		function output = induced_exogenous(A,source,target,nodelist,name,varargin)
			% 
			% Endogenous: standard measure of node centrality
			% Exogenous: indirect contribution to the centrality of other nodes
			
			if(~exist('normalized','var'))
				normalized = false;
			end
			
			if(~exist('name','var'))
				name = 'closeness';
			end
			
			if(~exist('metricfun'))
				metricfun = @(x,source,target)(communicability.invariant(x,source,target,name,normalized)); 							
			end
			
			[invariant end_centrality] = metricfun(A,source,target);
			
			centrality = sparse(zeros(length(A),1)); centrality_pos = centrality; centrality_neg = centrality;
			exogenous = sparse(zeros(length(A),1));
			endogenous = sparse(zeros(length(A),1));	
			induced = sparse(zeros(length(A),1));
			dyadic = spalloc(length(nodelist),length(A),1);
			C = [];
			for nn=1:length(nodelist)
				A_nn = A; 
				C(nn) = (length(A)-length(nodelist{nn}));
				A_nn(nodelist{nn},:) = 0; A_nn(:,nodelist{nn}) = 0;
				endogenous(nodelist{nn}) = end_centrality(nodelist{nn});
				[invariant_nn centrality_nn] = metricfun(A_nn,source,target);
				induced(nodelist{nn}) = invariant - invariant_nn;
				exogenous(nodelist{nn}) = (induced(nodelist{nn}) - sum(endogenous(nodelist{nn}))); 
				centrality_change = (end_centrality-centrality_nn)./end_centrality;
				centrality_change(nodelist{nn}) = 0; % 
				dyadic(nodelist{nn},:) = repmat(reshape(centrality_change, [1 length(A)]), [length(nodelist{nn}) 1]);							
				centrality(nodelist{nn}) = sum(abs(centrality_change));
				centrality_pos(nodelist{nn}) = sum(centrality_change.*(centrality_change>0));
				centrality_neg(nodelist{nn}) = sum(centrality_change.*(centrality_change<0));
			end
			
			output.normalization = C;	
			output.endogenous = endogenous;
			output.exogenous = exogenous;
			output.induced = induced;
			output.centrality = centrality;
			output.centrality_pos = centrality_pos;
			output.centrality_neg = centrality_neg;
			output.dyadic = dyadic;
			output.nodelist = nodelist;
			output.name = name;
			output.normalized = normalized;
			
		end
		
		function output = closeness(A,K,varargin)
			% K should be communicability.distance
			if(isempty(K))
				K = communicability.distance(A);
			end	
			p = size(K,1);
			
			switch nargin
			case 3
				subset = varargin{1};
				normalized=false;
			case 4
				subset = varargin{1};
				normalized=varargin{2};
			otherwise
				subset = 1:p;
				normalized=false;
			end
			
			% % standard closeness
			% farness = sum(K(:,subset),2)/p;			
			% closeness = 1./farness;
			% reverse closeness
			closeness = (length(subset) - sum(K(:,subset),2))/length(subset);
			
			if(normalized)
				output = closeness/max(closeness);
			else
				output = closeness;
			end
			
		end
		
		function output = global_centrality(varargin)
			
			warning('Not implemented')
			
		end
				
		function output = edge_metrics(A,varargin)
			
			if(nargin==2)
				method = varargin{1}; 
			else
				method = 'similarity'
			end
			
			if(isempty(G))
				switch method
				case 'similarity'
					G = communicability.communicability_matrix(A);
				case 'dissimilarity'
					G = communicability.distance(A);
				end
			end	
			
			G(find(eye(size(G,1)))) = 0;			
			output = G; 
				
		end
		
		
		function output = nodal_metrics(A,varargin)
			
			
		end
	
		function output = global_metrics(A,varargin)
			
			
		end
	
		function output = metrics(A,varargin)
			
			metrics = []; 
			labels = {}; 
			
			G = communicability.communicability_matrix(A); 
			K = communicability.distance(A,G); 
			
			metrics(:,1) = communicability.neighborhood(A,G,'subgraph'); 
			labels{1} = 'Com. Neighbor.';
			metrics(:,2) = communicability.degree(A,G); 
			labels{2} = 'Com. Degree';
			metrics(:,3) = communicability.closeness(A,K);
			labels{3} = 'Com. Closeness';
			
			if(nargin>=2)
				nodelist = varargin{1}; 				
			else
				nodelist = [1:size(A,1)]; 
			end
			
			metrics(:,4) = communicability.betweenness(A,nodelist);
			labels{4} = 'Com. Betweenness'; 
		
			output.metrics = metrics;
			output.labels = labels;
			output.type = {'Neighbor.', 'Degree', 'Closeness', 'Betweenness'}; 
		
			output.metrics = metrics;
			output.labels = labels;
			
		end
	
	
		function output = reduced_source_target(A,source,target)
		
			p = length(A);
			old_idx = 1:p;
			include_nodes = ~(sum(A)==0);
			new_idx = 1:sum(include_nodes);
			new2old = find(include_nodes);
			
			A_reduced = A(include_nodes,include_nodes);		
			[sA sB] = ismember(source,new2old);
			new_source = sB(sB~=0);
			[tA tB] = ismember(target,new2old);
			new_target = tB(tB~=0);
						
			output.A = A_reduced;
			output.node_idx = find(include_nodes);
			output.source = new_source;
			output.target = new_target;
			
		end
		
		function output = normalize_helper(centrality,method)
			
			if(~exist('method','var'))
				method = 'max';
			elseif(isempty(method))
				method = 'max'
			end
			
			Cpos = 0;
			Cneg = 0;
			centrality_pos = centrality.*(centrality>0); 
			centrality_neg = centrality.*(centrality<0);
			
			switch method
				
			case 'max'
				Cpos = norm(centrality_pos,'Inf');	
				Cneg = norm(-centrality_neg,'Inf');			
			case 'frob'
				C = norm(centrality,2);
				C_pos = C; C_neg = C;
			case 'abs'
				C = norm(centrality,1);
				C_pos = C; C_neg = C;				
			end
			
			if(C_pos==0)
				C_pos = 1;
			end
			if(C_neg==0)
				C_neg = 1;
			end
			output = centrality_pos/C_pos + centrality_neg/C_neg;
			
		end
	
	end
	
end