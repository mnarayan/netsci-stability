function output = standardize_centrality(centrality,method)
	
	if(~exist('method','var'))
		method = 'max';
	elseif(isempty(method))
		method = 'max';
	end
	
	C_pos = 0;
	C_neg = 0;
	centrality_pos = centrality.*(centrality>0); 
	centrality_neg = centrality.*(centrality<0);
	
	switch method
		
	case 'max'
		Cpos = norm(centrality_pos,'inf');	
		Cneg = norm(-centrality_neg,'inf');			
	case 'frob'
		C = norm(centrality,2);
		C_pos = C; C_neg = C;
	case 'abs'
		C = norm(centrality,1);
		C_pos = C; C_neg = C;	
	case 'zscore'
		centrality_pos = zscore(centrality_pos); 
		centrality_neg = zscore(centrality_neg); 
		C_pos = 1; 
		C_neg = 1;					
	end
	
	if(C_pos==0)
		C_pos = 1;
	end
	if(C_neg==0)
		C_neg = 1;
	end
	output = centrality_pos/C_pos + centrality_neg/C_neg;
	
end