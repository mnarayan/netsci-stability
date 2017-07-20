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
			D(D==0) = 1; 
			normalize = true;
			if(normalize)
				output = eye(p) - diag(1./sqrt(D))*A*diag(1./sqrt(D));
			else
				output = diag(D) - A; 
			end
			output = (output + output')/2;
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
		
		function output = neighborhood(A,G,varargin)
			% Eigenvector
			% Compare with walk probability matrix here
			% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4897067/
			
			if(isempty(G))
				G = currentflow.inverse_laplacian_matrix(A); 
			end
			p = size(G,1);
			% W = zeros(p,p);
			% W = A;
			% W(find(eye(p))) = 0.0;
			% D = sum(W,2);
			% D(D==0) = 1.0;
			% D_alt = sqrtm(diag(D));
			% W = D_alt + D_alt^(-1)*W*D_alt^(-1);
			%centrality = eigenvector_centrality_und(G);			
			%centrality = bonanich_power_centrality(W',struct('beta',.8));
			centrality = resolvent_centrality(G); 	
			output = centrality/sum(centrality); 
			
		end
		
		function [output] = conductance(A,Ci,varargin)
			% conductance of a cut
			%
			% Reference: 
			% Hierarchical Directed Spectral Graph Partitioning
			% by D. Gleich
			
			p = size(A,1); 
			A = triu(A,1); 
			output = zeros(p,1); 
			n_com = length(unique(Ci)); 
			
			for community_no = 1:n_com
				node_no = find(Ci==community_no); 
				other_nodes = setdiff([1:p],node_no); 
				
				vol_nodes = sum(sum(A(node_no,:),2));
				vol_othernodes = sum(sum(A(other_nodes,:),2)); 
				output(node_no) = sum(A(node_no,other_nodes))/min(vol_nodes,vol_othernodes); 
			end	
			
		end
		
		function [output varargout] = degree(A,G,varargin)
			% Degree and closeness are virtually identical 
			
			if(isempty(G))
				G = currentflow.inverse_laplacian_matrix(A); 
			end
			p = size(G,1); 
			isSigned = true;
			
			% K = triu(1./(G+eye(p)),1);
			% K = K + K';
			
			output = zeros(p,1); 	
			% rescale from 0 to 1
			% output = bsxfun(@minus,max(output),output)./(max(output)-min(output)); 
			% rescale by sum
			% output = output/sum(output); 
			
			if(isSigned)
				pos_centrality = sum(G.*(G>0)); 
				neg_centrality = abs(sum(G.*(G<0))); 
				output = pos_centrality/max(max(pos_centrality),1) + neg_centrality/max(max(neg_centrality),1); 
			else
				output = sum(G); 
			end
			
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
				for tt=1:ss
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
		
		
		function output = betweenness(A,G,varargin)
			% K should be currentflow.inverse_laplacian_matrix
			p = size(A,1); 
			output = zeros(p,1);

			if(exist('G','var'))
				if(isempty(G))
					G = currentflow.inverse_laplacian_matrix(A); 
				end
			else
				G = currentflow.inverse_laplacian_matrix(A); 
			end
			
			if(nargin==3)
				nodelist = varargin{1}; 
			else
				nodelist = 1:p;
			end 
			
			if(isempty(nodelist)|nodelist(1) == 0)
				disp('Skipping Betweenness')
			else
				for ii=[nodelist]
					%sprintf('Computing betweenness for node %d',ii)
					output(ii) = sum(sum(currentflow.flowmatrix(A,G,ii)))/nchoosek(p,2); 
				end
			end
			
            try
			    assert(sum(output<0)==0,'Negative values of Betweennes')
            catch me
                disp(me)
            end
			% disp('No. of negative values')
			% sum(output<0)
			
			% pos_output = output.*(output>0);
			% pos_output = pos_output/max(max(pos_output),1);
			% neg_output = abs(output.*(output<=0));
			% neg_output = neg_output/max(max(neg_output),1);
			%
			% output = pos_output + neg_output;
			
		end
		
		function output = closeness(A,K,varargin)
		% CURRENTFLOW.CLOSENESS 	
		% Inputs: 	
		% K should be currentflow.distance
		
            if(exist('K','var'))
    			if(isempty(K))
    				K = currentflow.distance(A); 
    			end	
            else
                K = currentflow.distance(A);
            end
			p = size(K,1);
			farness = sum(K)/p;
			closeness = 1./farness;
			closeness(closeness==Inf) = 0;
			output = closeness;
			
		end
		
		function output = global_efficiency(A,K,varargin)
		% CURRENTFLOW.GLOBAL_EFFICIENCY	
			% K should be currentflow.distance
			
			if(isempty(K))
				K = currentflow.distance(A); 
			end	
			p = size(K,1);
			distmatrix = K(triu(K,1)>0);
			if(length(distmatrix)<nchoosek(p,2))
				warning('efficiency does not have full upper triangular');
			end			
			efficiency = sum(sum(1./distmatrix))/nchoosek(p,2);
			output = efficiency;
		end
		
		function output = local_efficiency(A,varargin)
		% CURRENTFLOW.LOCAL_EFFICIENCY	
						
			p = size(A,1); 
			centrality = zeros(p,1);
			
			if(nargin>=2)
				nodelist = varargin{1}; 
			else
				nodelist = 1:p;
			end 
			
            
			if(isempty(nodelist)|nodelist{1} == 0)
				disp('Skipping Local. Efficiency')
			else
				for ii=1:length(nodelist)
					tmpA = zeros(p-1,p-1);
					idx = setdiff(1:p,nodelist{ii}); 
					tmpA = A(idx,idx);
					distmatrix = currentflow.distance(tmpA);
					centrality(nodelist{ii}) = ...
                         currentflow.global_efficiency(tmpA,distmatrix);
				end
			end
			output.centrality = centrality;
		end
		
		function output = metrics(A,varargin)
			
			std_method = 'max'
			metrics = []; 
			labels = {}; 
			
			G = currentflow.inverse_laplacian_matrix(A); 
			K = currentflow.distance(A,G); 
			
			metrics(:,1) = currentflow.neighborhood(A,G); 
			metrics(:,1) = standardize_centrality(metrics(:,1),std_method);
			labels{1} = 'CFlow Neighbor.';
			metrics(:,2) = currentflow.degree(A,K); 
			metrics(:,2) = standardize_centrality(metrics(:,2),std_method);
			labels{2} = 'CFlow Degree';
			metrics(:,3) = currentflow.closeness(A,K);
			metrics(:,3) = standardize_centrality(metrics(:,3),std_method);	 
			labels{3} = 'CFlow Closeness';
			
			if(nargin>=2)
				nodelist = varargin{1}; 				
			else
				nodelist = [1:size(A,1)]; 
			end
			metrics(:,4) = currentflow.betweenness(A,G,nodelist);
			metrics(:,4) = standardize_centrality(metrics(:,4),std_method);	
			labels{4} = 'CFlow Betweenness'; 
		
			if(nargin>=2)
				nodelist = varargin{1}; 				
			else
				nodelist = [1:size(A,1)]; 
			end
			metrics(:,5) = currentflow.local_efficiency(A,nodelist);
			metrics(:,5) = standardize_centrality(metrics(:,5),std_method);	
			labels{5} = 'CFlow LocEfficiency'; 
		
			output.metrics = metrics;
			output.labels = labels;
			output.type = {'Neighbor.',...
						   'Degree', ...
						   'Closeness', ...
						   'Betweenness', ...
						   'Loc Efficiency'}; 
			
		end
	
	end
	
end