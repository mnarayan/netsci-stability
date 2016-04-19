beta = .999;

% Test matrix
A = zeros(3,3);
A(1,3) = 1; A(2,3) = 1;
[hubs authorities eig_centrality] =  hits(A);
[h_pgrank] = pagerank_centrality(A, beta); [a_pgrank] = pagerank_centrality(A',beta);
centrality2centralization(eig_centrality,'eigenvector',1)

assert(~any(hubs - [.5 .5 0]'),'Hubs are incorrect for [0 0 1; 0 1 0; 0 0 0]');
assert(~any(authorities - [0 0 1]'),'Authorities are incorrect for [0 0 1; 0 1 0; 0 0 0]');
clear;


% Problematic Tests (Springer chapter)
% p1 <- p_2 -> p_3 -> p_4
A = zeros(4,4);
A(2,1) = 1; A(2,3) = 1; A(3,4) = 1;
[hubs authorities eig_centrality] =  hits(A);
[hubs authorities] =  rand_hits(A,1-beta);
[h_pgrank] = pagerank_centrality(A, beta); [a_pgrank] = pagerank_centrality(A',beta);
centrality2centralization(hubs,'hub',1)
centrality2centralization(authorities,'authority',1)
centrality2centralization(h_pgrank,'pagerank',1)
centrality2centralization(a_pgrank,'pagerank',1)
assert(~any(hubs - [0 1 0 0]'),'Hubs are incorrect');
assert(sum(authorities - [.5 0 .5 0]')< 1e-3,'Authorities are incorrect ');


% Problematic Tests (Springer chapter)
% p1 -> p_2; p1->p3; p2 -> p_3 -> p_4
A = zeros(4,4);
A(1,2) = 1; A(1,3) = 1;  A(2,3) = 1; A(3,4) = 1;
[hubs authorities eig_centrality] =  hits(A);
[hubs authorities] =  rand_hits(A,1-beta);
[h_pgrank] = pagerank_centrality(A, beta); [a_pgrank] = pagerank_centrality(A',beta);

% Cornell Example
A = zeros(4,4);
A(1,2) = 1; A(1,3) = 1;
A(3,4) = 1; A(3,2) = 1;
A(4,1) = 1;
[hubs authorities eig_centrality] =  hits(A);
[h_pgrank] = pagerank_centrality(A, beta); [a_pgrank] = pagerank_centrality(A',beta);
centrality2centralization(hubs,'hub',1)
centrality2centralization(authorities,'authority',1)
centrality2centralization(h_pgrank,'pagerank',1)
centrality2centralization(a_pgrank,'pagerank',1)
assert(~any(hubs - [.5 0 .5 0]'),'Hubs are incorrect');
assert(sum(authorities - [0 .5 .25 .25]')< 1e-3,'Authorities are incorrect ');
clear;


% Stanford NLP Example
A = zeros(7,7);
A(1,3) = 1; A(2,3) = 1; A(2,2) = 0;
A(3,1) = 1; A(3,3) = 0; A(3,4) = 2;
A(4,4) = 0; A(4,5) = 1; A(5,7) = 1;
A(6,6) = 0; A(6,7) = 1; A(7,4) = 2; A(7,5) = 1;, A(7,7) = 0;
[hubs authorities eig_centrality] =  hits(A);
[h_pgrank] = pagerank_centrality(A, beta); [a_pgrank] = pagerank_centrality(A',beta);
centrality2centralization(hubs,'hub',1)
centrality2centralization(authorities,'authority',1)
centrality2centralization(h_pgrank,'pagerank',1)
centrality2centralization(a_pgrank,'pagerank',1)
assert(sum([hubs(3) hubs(7)] - [0.4669 0.4742])< 1e-3,'Hubs are incorrect');
assert((authorities(4) - .6531)<1e-3,'Authorities are incorrect');
assert(sum(eig_centrality([3:4 7]) - [.1886 .1886 .1886]')< 1e-3,'Eig Centrality is incorrect');
clear;


% White & Smyth Undirected.
G = zeros(10,10)
G(1,2) = 1; G(2,3) = 1; G(3,1) = 1;
G(4,6) = 1; G(5,6) = 1; G(6,7) = 1;
G(1,4) = 1; G(4,5) = 1; G(5,10) = 1; G(10,3) = 1;
G(10,8) = 1; G(8,7) = 1; G(7,9) = 1; G(8,9) = 1;G(9,2) = 1;
G = G+G';
A = G;
[hubs authorities eig_centrality] =  hits(A);
[h_pgrank] = pagerank_centrality(A, beta); [a_pgrank] = pagerank_centrality(A',beta);
centrality2centralization(eig_centrality,'eigenvector',1);


% Directed/Asymmetric
G = zeros(10,10)
G(1,3) = 1; G(2,1) = 1; G(2,3) = 1; G(3,10) = 1; G(10,5) = 1; G(5,4) = 1; G(4,1) = 1; G(4,6) = 1; G(6,7) = 1; G(9,2) = 1;
G(5,6) = 1; %G(6,5) = 1;
G(8,10) = 1;%G(10,8) = 1;
G(7,8) = 1; %G(8,7) = 1;
G(9,8) = 1; %G(8,9) = 1;
G(9,7) = 1; %G(7,9) = 1;
A = G;
[hubs authorities eig_centrality] =  hits(A);
[h_pgrank] = pagerank_centrality(A, beta);
centrality2centralization(hubs,'hub',1)
centrality2centralization(authorities,'authority',1)
centrality2centralization(eig_centrality,'eigenvector',1)
centrality2centralization(h_pgrank,'pagerank',1)


%  White & Smyth, Directed.
G = zeros(10,10)
G(1,3) = 1; G(2,1) = 1; G(2,3) = 1; G(3,10) = 1;
G(10,5) = 1; G(5,4) = 1; G(4,1) = 1; G(4,6) = 1; G(5,6) = 1; G(6,5) = 1;G(6,7) = 1;
G(9,2) = 1;
G(10,8) = 1; G(8,10) = 1;
G(7,8) = 1; G(8,7) = 1;
G(8,9) = 1; G(9,8) = 1;
G(7,9) = 1; G(9,7) = 1;
A = G;
[hubs authorities eig_centrality] =  hits(A);
[h_pgrank] = pagerank_centrality(A, beta);
centrality2centralization(hubs,'hub',1)
centrality2centralization(authorities,'authority',1)
centrality2centralization(eig_centrality,'eigenvector',1)
centrality2centralization(h_pgrank,'pagerank',1)
