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
	
		function [output varargout] = total(A,varargin)
			% total communicability
			
			K = communicability.communicability_matrix(A); 
			K(find(eye(size(K,1)))) = 0;			
			output = sum(K); 
			
		end
	
		function output = subgraph_centrality(A,varargin)
			% subgraph centrality
			
			K = communicability.communicability_matrix(A); 			
			output = diag(K)'; 
			
		end
		
		function output = betweenness(varargin)
			
		end
		
		function output = closeness(varargin)
			
		end
		
		function output = global_centrality(varargin)
			
			
		end
		
	
		function output = metrics(A,varargin)
			
			
		end
	
	
	end
	
end