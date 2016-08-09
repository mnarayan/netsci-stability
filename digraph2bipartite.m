function A = digraph2bipartite(DG)
	
	
	assert(sum(sum(abs(DG-DG')))>1,'Adjacency matrix may not be assymetric');
	
	p = size(DG,1); 
	A = zeros(2*p,2*p);
	
	A(1:p,p+1:2*p) = DG;
	A(p+1:2*p,1:p) = DG';
	
	tmp_eig = eig(A);
	assert(~any(imag(tmp_eig(:))), 'Not symmetric');
	
	% Test bi-partite: This gives back hubs and authorities
	% B = digraph2bipartite(A);
	% [V D] = eig(B);
	% [sorted_D idx_und]  = sort(diag(abs(D(1:size(A,1),1:size(A,1)))),'descend');
	% if(visible)
	% 	disp('Largest two eigenvalues of symmetricized A ')
	% 	sorted_D(1:2)'
	% end
	% eig_centrality = V(size(A,1)+1:2*size(A,1),idx_und(1)); % ordered from highest to lowest
	% eig_centrality = eig_centrality /sum((eig_centrality));
	
end