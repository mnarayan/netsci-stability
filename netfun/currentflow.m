classdef currentflow
% All currentflow or random walk metrics
% Manjari Narayan
% Copyright 2017

	methods(Static)
	
		function output = inverse_laplacian_matrix(A,varargin)
			% inverse_laplacian_matrix.m
			% function output = inverse_laplacian_matrix(A,varargin)
			% 
			% inverse_laplacian(u,v) = 
			% 
			% INPUTS
			% A is an adjacency matrix for a weighted or unweighted graph
			% optional:
			%			% 
			% OUTPUTS:
			
			L = currentflow.laplacian(A); 
			output = pinv(L); 
		
		end

		function [output varargout] = laplacian(A, varargin)
		
			p = size(A,1); 
			A(find(eye(p))) = 0; 
				
			D = sum(abs(A)); 
			output = diag(D) - A; 
			
		end

		function [output varargout] = distance(A, varargin)
			% currentflow.distance (Resistance Distance)
			% function output = currentflow.distance(A,varargin)
			% 
			% resistance(u,v) = 
			% 
			% INPUTS
			% A is an adjacency matrix for a weighted or unweighted graph
			% optional:
			%		
			% OUTPUTS:
		
			if(nargin==1)
				K = currentflow.inverse_laplacian_matrix(A); 
			else
				K = varargin{1}; 
			end
			p = size(K,1); 			
			D = diag(K);
			
			output = zeros(size(K)); 
			for ii=1:p
				for jj=1:ii
					output(ii,jj) = D(ii) + D(jj) - 2*K(ii,jj); 
				end
			end

			output = output + output'; 
		
			% output.function = 'resistance_distance';
			% output.matrix = resistance;
			% output.input = laplacian_output;
		
		end
		
		function output = neighborhood(A,K, varargin)
			% Eigenvector | Katz | Bonanich
			
			
			
			% p = size(A,1);
			% A(find(eye(p))) = 0;
			% D = sum(abs(A));
			% D = reshape(D, [p 1]);
			% centrality = K*D;
			
			%centrality = eigenvector_centrality_und(A);
			centrality = bonanich_power_centrality(K); 	
			
			output = centrality; 
			
			
		end
		
		
		function [output varargout] = degree(A,varargin)
			% total communicability
			
			K = currentflow.inverse_laplacian_matrix(A); 
			K(find(eye(size(K,1)))) = 0;
			output = sum(K); 
			
			output = output/sum(output); 
			
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
				for tt=1:ss
					res = K(node,ss) - K(node,tt) + K(not_node,tt) - K(not_node,ss);
					output(ss,tt) = .5*sum(A(node,not_node).*res');
				end
			end
			
		end
		
		
		function output = betweenness(A,K,varargin)
			% K should be currentflow.distance
			p = size(A,1); 
			output = zeros(p,1);
			
			if(nargin==3)
				nodelist = varargin{1}; 
			else
				nodelist = 1:p;
			end 
			
			for ii=nodelist
				%sprintf('Computing betweenness for node %d',ii)
				output(ii) = sum(sum(currentflow.flowmatrix(A,K,ii)))/nchoosek(p,2); 
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
			% K should be currentflow.distance
			
			p = size(K,1);
			farness = sum(K)/p;
			closeness = 1./farness;
			output = closeness;
			
		end
		
		
		function output = metrics(A,varargin)
			
			metrics = []; 
			labels = {}; 
			
			G = currentflow.inverse_laplacian_matrix(A); 
			K = currentflow.distance(A,G); 
			
			metrics(:,1) = currentflow.neighborhood(A,G); 
			labels{1} = 'CFlow Neighbor.';
			metrics(:,2) = currentflow.degree(A); 
			labels{2} = 'CFlow Degree';
			metrics(:,3) = currentflow.closeness(A,K); 
			labels{3} = 'CFlow Closeness';
			
			if(nargin==2)
				nodelist = varargin{1}; 				
				metrics(:,4) = currentflow.betweenness(A,G,nodelist);
			else	
				metrics(:,4) = currentflow.betweenness(A,G);
			end
			labels{4} = 'CFlow Betweenness'; 
		
			output.metrics = metrics;
			output.labels = labels;
			output.type = {'Neighbor.', 'Degree', 'Closeness', 'Betweenness'}; 
			
		end
	
	end
	
end