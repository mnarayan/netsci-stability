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
			
			output = zeros(size(G)); 
			for ii=1:p
				for jj=1:ii
					output(ii,jj) = D(ii) + D(jj) - 2*G(ii,jj); 
				end
			end
			
			output = output + output'; 
			
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
		
		
		function output = betweenness(A,K,varargin)
			% K should be
			p = size(A,1); 
			output = zeros(p,1);
			
			if(nargin==3)
				nodelist = varargin{1}; 
			else
				nodelist = 1:p;
			end 
			
			if(isempty(nodelist)|nodelist(1) == 0)
				disp('Skipping betweenness')
			else
				for ii=nodelist
					%sprintf('Computing betweenness for node %d',ii)
					output(ii) = sum(sum(communicability.flowmatrix(A,K,ii)))/nchoosek(p,2); 
				end
			end
			
			assert(sum(output<0)==0,'Negative values of Betweennes')
			% disp('No. of negative values')
			% sum(output<0)
			
			output = output/max(max(output),1.0); 
			% pos_output = output.*(output>0);
			% pos_output = pos_output/max(max(pos_output),1);
			% neg_output = abs(output.*(output<=0));
			% neg_output = neg_output/max(max(neg_output),1);
			%
			% output = pos_output + neg_output;
			
		end
		
		function output = closeness(A,K,varargin)
			% K should be communicability.distance
			p = size(K,1);
			farness = sum(K)/p;
			closeness = 1./farness;
			output = closeness/max(closeness);
			
		end
		
		function output = global_centrality(varargin)
			
			
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
			
			metrics(:,4) = communicability.betweenness(A,G,nodelist);
			labels{4} = 'Com. Betweenness'; 
		
			output.metrics = metrics;
			output.labels = labels;
			output.type = {'Neighbor.', 'Degree', 'Closeness', 'Betweenness'}; 
		
			output.metrics = metrics;
			output.labels = labels;
			
		end
	
	
	end
	
end