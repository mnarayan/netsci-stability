function output = relative_importance(metric)
	% function output = relative_importance(metric)
	% 
	% Returns a table with ranks
	
	useMatlab = (exist('tiedrank')>=2);
	p = length(metric);
	metric = reshape(metric,[p 1]); 
	
	if(useMatlab)
		rank_metric = p+1 - tiedrank(metric,0); 
	end
	
	output = array2table(metric, 'VariableNames', {'rawmetric'});
	output.Properties.RowNames = cellfun(@num2str,num2cell(1:p),'UniformOutput',false);
	output.rank = rank_metric;
	
end
