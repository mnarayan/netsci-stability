classdef communicability
% All communicability related metrics
% Manjari Narayan
% Copyright 2017

	methods(Static)
	
		function [output varargout] = communicability_matrix(A,varargin)
			% communicability_matrix.m
			% function output = communicability(A,varargin)
			% 
			% communicability(u,v) = 
			% 
			% INPUTS
			% A is an adjacency matrix for a weighted or unweighted graph
			% optional:
			%			% 
			% OUTPUTS:
			output = expm(A); 
		
		end
	
		function [output varargout] = distance(A, varargin)
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
			K = communicability.communicability_matrix(A); 
			p = size(K,1); 			
			D = diag(K);
			
			output = zeros(size(K)); 
			for ii=1:p
				for jj=1:ii
					output(ii,jj) = D(ii) + D(jj) - 2*K(ii,jj); 
				end
			end
			
			output = output + output'; 
			
			
		end
	
		function [output varargout] = degree(A,K,varargin)
			% total communicability
			
			K = communicability.communicability_matrix(A); 
			K(find(eye(size(K,1)))) = 0;			
			output = sum(K); 
			
			output = output/max(output); 
			
		end
	
		function output = subgraph(A,K,varargin)
			% subgraph centrality
			
			output = diag(K)'; 
			
			output = output/sum(output);
			
		end
		
		function output = neighborhood(A,K, varargin)
			% Eigenvector
			
			
			
			centrality = eigenvector_centrality_und(K);			
			
			output = centrality; 
			
			
		end
		
		function output = flowmatrix(A,K,node,varargin)
			% F(s,t)i=12∑jAi,j|G+i,s − G+i,t + G+j,t − G+j,s|
			
			p = size(A,1); 
			output = zeros(p,p); 
			
			if(isempty(node)|nargin<3)
				error('node argument missing')
			end
			
			K(find(eye(size(K,1)))) = 0;		
			not_node = setdiff([1:p],node);
			for ss=1:p
				for tt=ss:p
					res = K(node,ss) - K(node,tt) + K(not_node,tt) - K(not_node,ss);
					output(ss,tt) = .5*sum(A(node,not_node).*res');
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
			
			for ii=nodelist
				%sprintf('Computing betweenness for node %d',ii)
				output(ii) = sum(sum(communicability.flowmatrix(A,K,ii)))/nchoosek(p,2); 
			end
			
			disp('No. of negative values')
			sum(output<0)
			
			pos_output = output.*(output>0); 
			pos_output = pos_output/max(pos_output); 
			neg_output = abs(output.*(output<=0)); 
			neg_output = neg_output/max(neg_output); 
			
			output = pos_output + neg_output;
			
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
		
	
		function output = metrics(A,varargin)
			
			metrics = []; 
			labels = {}; 
			
			G = communicability.communicability_matrix(A); 
			K = communicability.distance(A); 
			
			metrics(:,1) = communicability.subgraph(A,G); 
			labels{1} = 'Com. Neighbor.';
			metrics(:,2) = communicability.degree(A); 
			labels{2} = 'Com. Degree';
			metrics(:,3) = communicability.closeness(A,K);
			labels{3} = 'Com. Closeness';
			
			if(nargin==2)
				nodelist = varargin{1}; 				
				metrics(:,4) = communicability.betweenness(A,G,nodelist);
			else	
				metrics(:,4) = communicability.betweenness(A,G);
			end
			labels{4} = 'Com. Betweenness'; 
		
			output.metrics = metrics;
			output.labels = labels;
			output.type = {'Neighbor.', 'Degree', 'Closeness', 'Betweenness'}; 
			
			% metrics(:,3) = currentflow.neighborhood(A,K);
			% labels{3} = 'CurrentFlow_Neighborhood';
		
			output.metrics = metrics;
			output.labels = labels;
			
		end
	
	
	end
	
end