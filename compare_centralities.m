classdef compare_centralities

		methods(Static)

				function output = run(metrics,varargin)
	
					ranks = compare_centralities.rank_centrality(metrics); 
					plotmatrix_figh = compare_centralities.plotmatrix(ranks); 
					
					output.ranks = ranks;
					output.plotmatrix_figh = plotmatrix_figh;
					
					
				end

				function output = rank_centrality(metrics)
					
					n_metrics = size(metrics,2); 
					ranks = [];
										
					for metric_no = 1:n_metrics
						
						tmp_rank = relative_importance(round(metrics(:,metric_no),2)); 
						ranks = horzcat(ranks,table2array(tmp_rank(:,'rank'))); 
						
					end
					
					output = ranks;
					
				end

				function output = dvalues(ranks)
					
					% n_metrics = size(ranks,2);
					% output = zeros(n_metrics,n_metrics);
					%
					% for ii = 1:n_metrics
					% 	for jj = ii:n_metrics
					% 		output(ii,jj) =
					% 	end
					% end
					
					
				end


				function figh = plotmatrix(ranks,labels)
					
				figh = scatter(ranks(:,1),ranks(:,2),50);
				xlim([0 size(ranks,1)]); 
				ylim([0 size(ranks,1)]);  
					
				end
				
				function figh = plot_ranklines(ranks,labels)
	
	
				end
		end
end