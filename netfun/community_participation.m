function [degCi P] =community_participation(Sigma,Ci,varargin)
%COMMUNITY_PARTICIPATION   Community specific participation
%
%  metrics = community_participation(Sigma, Ci, opts);
%
%  This function computes the sum or weighted sum of all connections between each node and some set of nodes that belong to a community.
%  Within-Community case: If node in community A, the second set of nodes are all other nodes of community A
%  Between-Community case: If node in community A, the second set of nodes are all nodes in a community not A
%  Participation Coefficient: Nodes allegience to a community as measured by  [Between-Community degree]/[Total degree]
%
%   Inputs:     Sigma,      A sample correlation matrix (binary/weighted | saturated/thresholded)
%               Ci,     	A vector of integer values,  with integer indicating community membership
% 							opts,
%
%   Output:     metrics
%       				metrics.within		nodes x 1 where each row measure degree wrt to community of respective node
% 							metrics.between  	nodes x max(Ci)-1 where each row measures degree wrt to all other communities.
% 							metrics.participation 	
%
% 
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900 & Brain Connectivity Toolbox by Sporns & Rubinov
%
% SEE ALSO participation_coef, module_degree_zscore, degree_centrality, degree_und

Ci = 1*(Ci);
p = length(Sigma);                   %number of nodes
Sigma(find(eye(p))) = 0;
deg	= sum(abs(Sigma),2);             %degree
NeighborCi	= (Sigma~=0)*diag(Ci);   %neighbor community affiliation
degCi = zeros(p,1);


% Ensure only 1 community
Ci = (Ci~=0)*1;
degCi = sum(abs(Sigma).*(NeighborCi==1),2)./deg;
degCi(deg==0) = 0;
P = ones(p,1)- degCi;
