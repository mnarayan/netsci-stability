% Plot metric at a single node
figure; hold on;
for time=1:length(obs)
	subplot(length(obs),1,time);
	tmp_x = reshape(taus,[n_thresh 1]);
	tmp_y = squeeze(mean(netstat_nodes(:,time,node,:),1)); tmp_y = reshape(tmp_y,[n_thresh 1]);
	%plot(tmp_x,tmp_y,'linewidth',3); hold on;
	errorbar(tmp_x,tmp_y,squeeze(var(netstat_nodes(:,time,node,:),1))); hold on;
	title(sprintf('t = %3.1f',obs(time)),'fontsize',16);
	if(time==length(obs))
		xlabel('Varying threshold','fontsize',20');
	end
	if(time==ceil(length(obs)/2))
		ylabel('Network Metric','fontsize',20)
	end
	ylim([-75 75]);
	xlim([.08 .22]);
end
filename = ['Data/Discontinuity_Example2' '_' graphtype '_funnum' num2str(bct_num) '_Weighted' num2str(isWeighted)]
savefig(filename)

if(usePlotly)
	response = fig2plotly(gcf, 'offline',true, 'filename',filename,'fileopt', 'overwrite','strip','false');
	plotly_url = response.url;
	saveplotlyfig(response,filename);
end

fig1 = figure;
h = plot(tmp_x,squeeze(netstat_nodes(1,t_idx,node,:))'); hold on;
set(h,'LineWidth',3);
h= plot(tmp_x,squeeze(netstat_nodes(2,t_idx,node,:))');
set(h,'LineWidth',3);
h=plot(tmp_x,squeeze(netstat_nodes(3,t_idx,node,:))');
set(h,'LineWidth',3);
set(h,'LineWidth',3);
legend({['Monte Carlo Run 1, t=' num2str(obs(t_idx))], ...
				['Monte Carlo Run 2, t=' num2str(obs(t_idx))], ...
				['Monte Carlo Run 3, t=' num2str(obs(t_idx))]},'fontsize',20)
xlabel('Varying threshold','fontsize',20);
ylabel('Network Metric','fontsize',20)
switch bct_num
	case 1
	title('(Dis)Continuity of Betweenness Centrality for a Single Node','fontsize',24);
	case  2
	title('(Dis)Continuity of Clustering Coefficient for a Single Node','fontsize',24);
	otherwise
	title(['(Dis)Continuity of Metric: ' func2str(bct_funs{bct_num}) ', for a Single Node'],'fontsize',24);

end
filename1 = ['Data/Discontinuity_Example1' '_' graphtype '_funnum' num2str(bct_num) '_Weighted' num2str(isWeighted)]
savefig(filename1)

if(usePlotly)
	response = fig2plotly(gcf,'offline',true, 'filename',filename1,'fileopt', 'overwrite','strip','false');
	plotly_url = response.url;
	saveplotlyfig(response,filename1);
end
