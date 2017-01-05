classdef currentflow
% All currentflow or random walk metrics
% Manjari Narayan
% Copyright 2017

	methods(Static)
	
		function output = inverse_laplacian_matrix(A,varargin)
			% inverse_laplacian_matrix.m
			% function output = communicability(A,varargin)
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
		
			K = currentflow.inverse_laplacian_matrix(A); 
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
		
		function output = neighborhood(A,varargin)
			
			centrality = bonanich_power_centrality(A); 	
			output = centrality; 
			
		end
		
		
		function output = betweenness(varargin)
			
			
		end
		
		
		function output = closeness(varargin)
			
			
		end
		
		
		function output = metrics(A,varargin)
			
			
		end
	
	end
	
end