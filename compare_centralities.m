classdef compare_centralities

		methods(Static)

				function output = run(A,varargin)
	
					p = size(A,1); 
					if(nargin==2)
						nodelist = varargin{1}; 
					else
						nodelist = 1:p;
					end
				
					plotfigs = 0;
					metrics = []; 
					labels = {}; 
					centrality_type = {};
					centrality_labels = {'Neighborhood', 'Degree', 'Closeness','Betweenness'}; 
					traversal_labels = {'Geo.', 'Cflow.', 'Comm.'};
					addpath(genpath('../packages/BCT'));

					metrics(:,1) = eigenvector_centrality_und(A); 
					centrality_type{1} = @eigenvector_centrality_und;
					labels{1} = 'Eigenvector'; 
					
					metrics(:,2) = degrees_und(A); 
					labels{2} = 'Degree'; 
					centrality_type{2} = @degree_und;
					
					metrics(:,3) = efficiency_wei(A,1); 
					centrality_type{3} = @efficiency_wei; 
					labels{3} = 'Geo Closeness'; 
					
					metrics(:,4) = betweenness_wei(A); 
					centrality_type{4} = @betweenness_wei; 
					labels{4} = 'Geo Betweenness'; 

					metrics(1:2,:)
					tmp_output = currentflow.metrics(A,nodelist); 
					tmp_output.metrics(1:2,:)
					metrics = horzcat(metrics,tmp_output.metrics); 
					labels = horzcat(labels,tmp_output.labels);
					clear tmp_output
					
					tmp_output = communicability.metrics(A,nodelist); 
					tmp_output.metrics(1:2,:)
					metrics = horzcat(metrics,tmp_output.metrics); 
					labels = horzcat(labels,tmp_output.labels);
					clear tmp_output
					
					
					metric_idx = 1:size(metrics,2); 
					% metric_idx = [3,]
					
					ranks = compare_centralities.rank_centrality(metrics); 

					discriminability = compare_centralities.discriminability(metrics); 

					
					if(plotfigs)
						plotmatrix_figh = compare_centralities.plotmatrix(ranks,labels); 
					else
						plotmatrix_figh = [];
					end
					

					output.centrality_labels = centrality_labels;
					output.traversal_labels = traversal_labels;
					output.dims_idx = [1 2 3 4];
					output.neighborhood_idx = [1 5 9];  
					output.degree_idx = [2 6 10]; 
					output.closeness_idx = [3 7 11]; 
					output.betweenness_idx = [4 8 12]; 
					output.chosen = [output.degree_idx ...  
														output.neighborhood_idx ...
														output.closeness_idx ...
														output.betweenness_idx]; 
					output.metrics = metrics;
					output.labels = labels;
					output.ranks = ranks;
					
					output.plotmatrix_figh = plotmatrix_figh;
					
				end

				function output = discriminability(metrics)
					
					p = size(metrics,2); 
					output = zeros(1,p); 
					for ii=1:p
							output(ii) = centrality2centralization(metrics(:,ii),[],0); 
					end
					
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
					
					[h,AX,BigAx,P,PAx] = plotmatrix(ranks);
					% xlim([0 size(ranks,1)]);
					% ylim([0 size(ranks,1)]); 
					set(h,'MarkerSize',30);
					%set(AX,'xlim',[0 size(ranks,1)],'ylim',[0 size(ranks,1)]);
					
					set(gcf,'Position',[400 100 850 900])
					
					for ii=1:size(ranks,2)
						set(get(AX(end,ii),'XLabel'),'String',labels{ii}); 
						set(P(1,ii),'Visible','off')
					end
					
					for ii=1:size(ranks,2)
						set(get(AX(ii,1),'YLabel'),'String',labels{ii}); 
					end
					
					figh.h = h; 
					figh.ax = AX; 
					figh.bigax = BigAx;
					figh.p = P; 
					figh.pax = PAx; 
					
				end
				
				function figh = barplot(X,labels,varargin)
					
						addpath(genpath('../packages/gramm'))
					
						if(nargin==3)
							colors = varargin{1}; 
						else
							colors = [];
						end
						
						if(nargin==4)
							groups = varargin{1}; 
						else
							groups = [];
						end
					
						xdims = size(X);
						if(xdims(1)>1 && xdims(2)>1)
							error('Data matrix not support. Provide a vector')
						end 
					
						X = reshape(X, [max(xdims) 1]); 
						labels = reshape(labels, [max(xdims) 1]); 
					
						if(isempty(colors))
							g=gramm('x',labels,'y', X, 'color', labels);
						else
							g=gramm('x',labels,'y', X, 'color', colors);
						end
						%%% 
						% Subdivide the data in subplots horizontally by region of origin using
						% facet_grid()
						if(mod(max(xdims),6)==0)
							cen_type = {'Closeness', 'Neighborhood'};
							grid_labels = reshape(repmat(cen_type,[3 1]), [max(xdims) 1]);
							g.facet_grid([],grid_labels);
							%g.facet_wrap(grid_labels,'ncols',2)
						
						end
						%%%
						g.geom_bar();
						%%%
						% % Plot linear fits of the data with associated confidence intervals
						% g.stat_glm();
						%%%
						% Set appropriate names for legends
						g.set_names('x', 'Type','y','Centralization (Mean/Max)','size',24);
						g.set_text_options('base_size',18,'label_scaling',1.4,'title_scaling',1.4,'big_title_scaling',2)
						%%%
						% Set figure title
						g.set_title('Discriminability of Centralities');
						g.no_legend()
						g.set_color_options('map','brewer_dark','lightness',30,'chroma',50); 
						g.set_order_options('x',0)
						if(max(X)<1)
							g.axe_property('ylim',[0 1]); 
						end
					
						%%%
						% Do the actual drawing
						figure('Position',[100 100 1000 600]);
						g.draw();
					
					
					
				end
				
				
				function figh = plot_ranklines(ranks,labels)
					
					
	
				end
		end
end