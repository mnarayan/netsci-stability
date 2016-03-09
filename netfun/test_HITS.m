% Test matrix
A = zeros(3,3);
A(1,3) = 1; A(2,3) = 1;
[hubs authorities eig_centrality] =  HITS(A);
assert(~any(hubs - [.5 .5 0]'),'Hubs are incorrect for [0 0 1; 0 1 0; 0 0 0]');
assert(~any(authorities - [0 0 1]'),'Authorities are incorrect for [0 0 1; 0 1 0; 0 0 0]');
clear;


A = zeros(4,4);
A(1,2) = 1; A(1,3) = 1;
A(3,4) = 1; A(3,2) = 1;
A(4,1) = 1;
[hubs authorities eig_centrality] =  HITS(A);
assert(~any(hubs - [.5 0 .5 0]'),'Hubs are incorrect');
assert(sum(authorities - [0 .5 .25 .25]')< 1e-3,'Authorities are incorrect ');
clear;


% Stanford NLP Example
A = zeros(7,7); 
A(1,3) = 1; A(2,3) = 1; A(2,2) = 0; 
A(3,1) = 1; A(3,3) = 0; A(3,4) = 2; 
A(4,4) = 0; A(4,5) = 1; A(5,7) = 1;
A(6,6) = 0; A(6,7) = 1; A(7,4) = 2; A(7,5) = 1;, A(7,7) = 0;
[hubs authorities eig_centrality] =  HITS(A);
assert(sum([hubs(3) hubs(7)] - [0.4669 0.4742])< 1e-3,'Hubs are incorrect');
assert((authorities(4) - .6531)<1e-3,'Authorities are incorrect');
assert(sum(eig_centrality([3:4 7]) - [.1886 .1886 .1886]')< 1e-3,'Eig Centrality is incorrect');
clear;