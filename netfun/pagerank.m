% pagerank.m
% function centrality = pagerank(A,beta)
% M = (beta*U + (1-beta)*A_row)'
%   where A_row is a row stochastic adjacency matrix
%   where U is the all ones matrix normalized by 1/(no. of nodes)
% Thus eigenvector of M is the pagerank centrality
%
% See "Stable Algorithms for Link Analysis" by Ng, Zheng and Jordan
% Manjari Narayan
% Copyright 2016
