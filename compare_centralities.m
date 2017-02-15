classdef compare_centralities

	methods(Static)

		function output = run(A,varargin)

			addpath(genpath('../packages/BCT'));
			plotfigs = 0;
			metrics = []; 
			labels = {}; 
			centrality_type = {};
			centrality_labels = {'Neighborhood', 'Degree', 'Closeness','Betweenness'}; 
			traversal_labels = {'Geo.', 'Cflow.', 'Comm.'};
			p = size(A,1); 
			
			if(nargin>=2)
				nodelist = varargin{1}; 
			else
				nodelist = 1:p;
			end
			if(nargin==3)
				traversal_list = varargin{2}; 
			else
				traversal_list = traversal_labels;
			end
		
			for traversal_no = 1:length(traversal_labels);

				switch traversal_no
					
				case 1
					disp('Geodesic Type')
					if(ismember('Geo.',traversal_list))
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
						labels = reshape(labels,[1 length(labels)]); 
					else
						metrics = zeros(p,length(centrality_labels)); 
						labels = {'Eigenvector', 'Degree', 'Geo. Closeness', 'Geo. Betweenness'}; 
						labels = reshape(labels,[1 length(centrality_labels)]); 
						
					end
				case 2	
					disp('Current Flow Type')		
					if(ismember('Cflow.',traversal_list))		
						tmp_output = currentflow.metrics(A,nodelist); 
						tmp_output.metrics(1:2,:)
						
					else
						tmp_output.metrics = zeros(p,length(centrality_labels)); 
						tmp_output.labels = repmat({''},[1 length(centrality_labels)]);
					end
					metrics = horzcat(metrics,tmp_output.metrics); 
					labels = horzcat(labels,tmp_output.labels);
					clear tmp_output
				case 3		
					disp(['Communicability Type'])					
					if(ismember('Comm.',traversal_list))		
						tmp_output = communicability.metrics(A,nodelist); 
						tmp_output.metrics(1:2,:)
					else
						tmp_output.metrics = zeros(p,length(centrality_labels)); 
						tmp_output.labels = repmat({''},[1 length(centrality_labels)]);
					end
					metrics = horzcat(metrics,tmp_output.metrics); 
					labels = horzcat(labels,tmp_output.labels);
					clear tmp_output
				end
			end
			
			metric_idx = 1:size(metrics,2); 
			% metric_idx = [3,]
			discriminability = compare_centralities.discriminability(metrics); 
			ranks = compare_centralities.rank_centrality(metrics); 
			dvalues = compare_centralities.dvalues(ranks); 
			[rank_agreement total_agreement] = compare_centralities.intersection_distance(ranks);
			%rank_pairs = compare_centralities.intersection_distance(ranks);  
			global_agreement = compare_centralities.global_agreement(rank_agreement); 
			
			if(plotfigs)
				plotmatrix_figh = compare_centralities.compare_centrality_types(dvalues,...
																	centrality_labels, ...
																    traversal_labels);
			else
				plotmatrix_figh = [];
			end
			

			output.centrality_labels = centrality_labels;
			output.traversal_labels = traversal_labels;
			output.traversal_groups = reshape(repmat(traversal_labels,[length(centrality_labels) 1]),[length(labels) 1]);
			output.centrality_groups = reshape(repmat(centrality_labels,[length(traversal_labels) 1]),[length(labels) 1]);
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
			output.dvalues = dvalues;
			output.rank_agreement = rank_agreement;
			output.node_agreement = total_agreement;
			output.global_agreement = global_agreement;
			output.discriminability = discriminability; 
			output.plotmatrix_figh = plotmatrix_figh;
			
		end


		function output = run_metric(A,metricfun,label,nodelist,varargin)
			% run_metric
			% function output = run_metric(A,metricfun,label,varargin)
			% 
			% INPUTS
			% 	- A 	 		A p times p adjacency matrix 
			% 	- metricfun 	A function handle to metric
			% 	- label			A string containing labelname of the metric
			% OUTPUTS
			% 
			% 
			
			narginchk(1,6);
			
			if(nargin==6)
				source = varargin{1}; 
				target = varargin{2}; 
			elseif(nargin>4)
				error('Input arguments not supported')
			end
			
			if(exist('metricfun','var'))
				fprintf('Calculating centrality statistics using function:%s \n ',func2str(metricfun));
			else
				metricfun = @resolvent_centrality;
				fprintf('Using default metric %s', metricfun)
			end
			
			if(~exist('label','var'))
				fprintf('No label provided. Using name:%s \n ',func2str(metricfun));
				label = func2str(metricfun);
			end
			
			p = size(A,1);
			if(~exist('source','var'))
				source = 1:p;
			end
			if(~exist('target','var'))
				target = 1:p;
			end
			
			output = compare_centralities.create_metric_struct();			
			fun_results = metricfun(A,nodelist,source,target); 
			metrics = fun_results.centrality;
			vec_metrics = full(metrics(metrics~=0)); vec_metrics = reshape(vec_metrics,[length(vec_metrics) 1]);
			size(vec_metrics)
			metrics
			discriminability = compare_centralities.discriminability(vec_metrics); 
			ranks(metrics~=0) = sparse(compare_centralities.rank_centrality(vec_metrics));
			
			metrics(nodelist{1})
			output.fun = metricfun;
			output.label = label;
			output.source = source;
			output.target = target;
			output.nodelist = nodelist;
			output.fun_results = fun_results;
			output.metrics = metrics;
			output.ranks = ranks;
			output.discriminability = discriminability;
			% output = setfield(output,{1},...
			% 				'fun',metricfun,...
			% 				'label',label,...
			% 				'source',source, ...
			% 				'target',target, ...
			% 				'nodelist',nodelist, ...
			% 				'fun_results',fun_results,...
			% 				'metrics',metrics,...
			% 				'ranks',ranks,...
			% 				'discriminability',discriminability);
			%			
			
		end

		function output = create_metric_struct(varargin)
		
			output = struct();
			output.fun = @resolvent_centrality; 
			output.fun_results = {};
			output.label = 'resolvent'; 
			output.metrics = [];
			output.discriminability = [];
			output.rank_pairs = [];
			output.source = [];
			output.target = [];
			output.nodelist = [];
		
		end

		function output = discriminability(metrics)
			% 
			% Compute discriminability of a network  using centralization of each centrality metric.
			% Centralization is high if some nodes have high centrality and others have a low centrality. 
			% 
			% INPUT
			% - metrics consists of n_nodes x n_metrics
			% OUTPUT
			% - output is a vector of length 1 x n_metrics
			% 
			
			p = size(metrics,2); 
			output = zeros(1,p); 
			for ii=1:p
					output(ii) = centrality2centralization(metrics(:,ii),[],0); 
			end
			
		end

		function output = rank_centrality(metrics)
			% 
			% Uses matlab's tiedrank function to rank nodes according to their centrality
			% 
			% INPUT
			% - metrics consists of n_nodes x n_metrics
			% OUTPUT
			% - output consists of n_nodes x n_metrics
			
			n_metrics = size(metrics,2); 
			ranks = [];
								
			for metric_no = 1:n_metrics
					
				  tmp_rank = relative_importance(metrics(:,metric_no)); 
				  ranks = horzcat(ranks,table2array(tmp_rank(:,'rank'))); 
				
			end
			
			output = ranks;
			
		end


		function output = topk(ranks,varargin)
			% Transform ranks to list of nodes within top-k
			% 
			% INPUT
			% - ranks consists of n_nodes x n_metrics
			% - k (optional): a number from 1 to n_nodes but is n_nodes by default
			
			if(nargin>=2)
				k = varargin{1}; 
			else
				k = 1:size(ranks,1);
			end
			
			output = zeros(size(ranks,1),size(ranks,2),length(k));
			
			for kk=1:length(k)
				output(:,:,k(kk)) = ranks<=k(kk);
			end
			
		end
		
		
		function [output varargout]= intersection_distance(ranks,varargin)
			% The intersection distance compare agreement or similarity between two different metrics within each top-k radius
			% 
			
			topk_matrix = compare_centralities.topk(ranks); 
			[p nmetrics k] = size(topk_matrix); 					

			output = zeros(p,nmetrics,nmetrics,k); 					
			for ii=1:size(ranks,2); 
				for jj=ii:size(ranks,2); 
						output(:,ii,jj,:) = squeeze(topk_matrix(:,ii,:).*topk_matrix(:,jj,:));
				end
			end
			
			kmatrix = repmat(fliplr([1:k]')/sum(1:k), [ 1 p nmetrics nmetrics]); 
			kmatrix = permute(kmatrix,[2 3 4 1]); 
			nodal_agreement = sum(output.*kmatrix,4);	
			node_idx = find(sum(sum(nodal_agreement>eps,2),3)~=0);
			
			node_struct = {};
			nargoutchk(1,2);
			if(nargout>=2)
				for ii=1:length(node_idx)
						node_struct{ii}.node_no = node_idx(ii); 
						node_struct{ii}.agreement = sparse(squeeze(nodal_agreement(node_idx(ii),:,:)));
				end			
				varargout{1} = node_struct;
			end
			
		end
		
		function output = global_agreement(agreement,varargin)
			% Agreement between centrality metrics across whole network
			[p nmetrics nmetrics k] = size(agreement);
			
			kmatrix = repmat(fliplr([1:k]')/sum(1:k), [ 1 p nmetrics nmetrics]); 
			kmatrix = permute(kmatrix,[2 3 4 1]); 
			norm_agreement = sum(agreement.*kmatrix,4);
			%norm_agreement = squeeze(agreement(:,:,:,end))/max(1,k);
			norm_constant = squeeze(sum(norm_agreement~=0,1)); 
			norm_constant(norm_constant==0) = 1.0;
			output = squeeze(sum(norm_agreement,1))./norm_constant; 
			
		end
		
		function output = dvalues(ranks)
			% 
			% Given any two samples (x,y) a d-value or stress-strength reliability capture P(x>y) as an estimate of population effect size that is independent of sample size.
			% This function seeks to offer a similar notion of effect size for centrality metrics or their rank equivalents. 
			% 
			% 
			% TODO: Re-implement d-values using the bootstrap
			
			if(length(size(ranks))==3)
				[p nmetrics B] = size(ranks); 
			else
				[p nmetrics] = size(ranks); 
				B = 1;
			end
			
			output = zeros(p,nmetrics);
			for bb=1:B 
				for ii=1:p
					node_ii = setdiff(1:p,ii);
					auc = zeros(1,nmetrics); 
					auc = mean(bsxfun(@le,ranks(ii,:,bb),squeeze(ranks(node_ii,:,bb))),1);
					output(ii,:,bb) = auc;
				end
			end
		end
		


		function figh = compare_centrality_types(dvalues,centrality_labels,traversal_labels,varargin)
			
			switch nargin
				
			case 4
				savefigs = varargin{1}; 
				filename = ['Example_' datestr(date)];
			case 5
				savefigs = varargin{1}; 
				filename = varargin{2}; 
			otherwise
				savefigs = 0; 
				filename = ''; 
			end
			
			pdf_list = {};
			figh = {};
			for ii=1:length(traversal_labels);
				dims_idx = 4*(ii-1) + [1:length(centrality_labels)]; 
				for jj=1:ii
					dims_idx2 = 4*(jj-1) + [1:length(centrality_labels)]; 
					if(ii==jj)
						figh{ii} = compare_centralities.plotmatrix(dvalues(:,dims_idx),...
																	[],...
																	centrality_labels,...
																	'Node AUC',...
																	traversal_labels{ii});
					else
						figh{ii} = compare_centralities.plotmatrix(dvalues(:,dims_idx),...
																	dvalues(:,dims_idx2),...
																	centrality_labels,...
																	'Node AUC',...
																	{traversal_labels{ii}, traversal_labels{jj}});
					end
					if(savefigs)
						pause(2);
						pdf_list = horzcat(pdf_list,['Data/' filename '_' num2str(ii) '_' num2str(jj) '.tif']);
						export_fig(['Data/' filename '_' num2str(ii) '_' num2str(jj) '.png']);
						
					end
				end
			end
			if(savefigs)
				%append_pdfs(['Data/' filename],pdf_list{:});
				%unix(['montage ' 'Data/' filename '_*.png' ' -tile 3x3' ' Data/' filename '.png']);
				unix(['convert ' 'Data/' filename '*.png'  ' Data/' filename '.pdf']);
				%unix(['rm Data/' filename '*.tif']); 
			end
	  end
		
		function figh = plotmatrix(X,Y,labels,varargin)
			
			addpath(genpath('../packages/gramm'))
			addpath(genpath('../packages/suplabel'));
		
			xdims = size(X);
			if(xdims(1)>1 & xdims(2)>1)
				warning('Data matrix not support. Provide a vector')
				X = reshape(X, [xdims(1)*xdims(2) 1]);
				labels = reshape(labels,[1 xdims(2)]);   
				labels = reshape(repmat(labels,[xdims(1) 1]), [xdims(1)*xdims(2) 1] ); 
			end
			
			if(nargin>=4 & ~isempty(varargin{1}))
				datatype = varargin{1}; 
			else
				datatype = 'Node AUC';
			end
			if(nargin>=5)
				traversaltypes = varargin{2}; 
			else
				traversaltypes = 'Traversal: ';
			end
			if(length(traversaltypes)==2)
				traversaltype1 = [traversaltypes{1} ' (X) ']; 
				traversaltype2 = [traversaltypes{2} ' (Y) ']; 
			end
			
			xdims = size(X);
			X = reshape(X, [max(xdims) 1]); 
			labels = reshape(labels, [max(xdims) 1]); 
			
			if(~isempty(Y))
				ydims = size(Y);
				if(ydims(1)>1 & ydims(2)>1)
					Y = reshape(Y, [ydims(1)*ydims(2) 1]);
				end
				ydims = size(Y); 
				Y = reshape(Y, [max(ydims) 1]);
			end 
			
			metric_names = unique(labels); 
			n_metrics = length(metric_names); 
			g = gramm();
			for ii=1:n_metrics
				for jj=ii:n_metrics
					if(ii==jj)
						ii_idx = strcmp(labels,metric_names{ii}); 
						if(isempty(Y))
							g(ii,ii) = gramm('x',X(ii_idx));
							g(ii,ii).stat_density('kernel','epanechnikov');
							g(ii,ii).set_title(metric_names{ii});
							g(ii,ii).set_names('x',datatype,'size',12);
							g(ii,ii).set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
						else
							g(ii,ii) = gramm('x',X(ii_idx),'y',Y(ii_idx));
							g(ii,ii).geom_point();
							g(ii,ii).stat_smooth();
							g(ii,ii).set_title(metric_names{ii});
							g(ii,ii).set_names('x', [ metric_names{ii} ' (X)'] , 'y', [ metric_names{ii} ' (Y)'],'size',12);
							g(ii,ii).set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
						end
						
					else
						ii_idx = strcmp(labels,metric_names{ii}); 
						jj_idx = strcmp(labels,metric_names{jj}); 
						
						if(isempty(Y))
							g(ii,jj) = gramm('x',X(ii_idx),'y',X(jj_idx));
							g(ii,jj).geom_point();
							g(ii,jj).stat_smooth();
							g(ii,jj).set_names('x', metric_names{ii}, 'y', metric_names{jj},'size',12);
							g(ii,jj).set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
						else
							g(ii,jj) = gramm('x',X(ii_idx),'y',Y(jj_idx));
							g(ii,jj).geom_point();
							g(ii,jj).stat_smooth();
							g(ii,jj).set_names('x', [ metric_names{ii} ' (X)'] , 'y', [ metric_names{jj} ' (Y)'],'size',12);
							g(ii,jj).set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
						end
						
						if(isempty(Y))
							g(jj,ii) = gramm('x',X(jj_idx),'y',X(ii_idx));
							g(jj,ii).geom_point();
							g(jj,ii).stat_smooth();
							g(jj,ii).set_names('x', metric_names{jj} , 'y', metric_names{ii},'size',12);
							g(jj,ii).set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
						else
							g(jj,ii) = gramm('x',X(jj_idx),'y',Y(ii_idx));
							g(jj,ii).geom_point();
							g(jj,ii).stat_smooth();
							g(jj,ii).set_names('x', [metric_names{jj} ' (X)'] , 'y', [ metric_names{ii} ' (Y)'],'size',12);
							g(jj,ii).set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
						end
					end
				end
			end
			
			if(isempty(Y))
				g.set_title(traversaltypes); 
			else
				g.set_title([traversaltype1 ' vs. ' traversaltype2]);
			end
			g.set_text_options('base_size',12,'label_scaling',1.35,'title_scaling',1.4,'big_title_scaling',1.8);
			
			figure('Position',[150 75 800 750]);
			g.draw();
			figh = g;
			
			if(isempty(Y))
				for ii=1:n_metrics
					set(get(figh(ii,ii).facet_axes_handles,'YLabel'),'String','pdf')
				end
			end
			
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
		
	end
end