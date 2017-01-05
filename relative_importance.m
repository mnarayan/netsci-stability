function output = relative_importance(metric)
	% function output = relative_importance(metric)
	% 
	% Returns a table with ranks
	
	p = length(metric);
	metric = reshape(metric,[p 1]); 
	
	rank_metric = tiedrank(metric,0); 
	
	output = array2table(metric, 'VariableNames', {'rawmetric'});
	output.Properties.RowNames = cellfun(@num2str,num2cell(1:p),'UniformOutput',false);
	output.rank = rank_metric;
	
end