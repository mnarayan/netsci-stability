% callnetworkx
% Calls the networkx current flow betweenness function;
% callnetworkx(A, isWeighted)
% callnetworkx(A, isWeighted, isScaled);

function centrality = callnetworkx(A, isWeighted,varargin);

isScaled = 0;
if(nargin>=3)
  isScaled = varargin{1};
end

try
	tbl = array2table(round(A,3));
catch
	tbl = array2table(A);
end
tbl.Properties.RowNames = tbl.Properties.VariableNames;
writetable(tbl,'tmp_A.txt','Delimiter','\t','WriteVariableNames',1,'WriteRowNames',1);

if(isWeighted)
  if(isScaled)
    unix(['python python/node.py -w -n tmp_A.txt tmp_output_A.txt']);
  else
    unix(['python python/node.py -w tmp_A.txt tmp_output_A.txt']);
  end
else
  if(isScaled)
    unix(['python python/node.py -n tmp_A.txt tmp_output_A.txt']);
  else
    unix(['python python/node.py  tmp_A.txt tmp_output_A.txt']);
  end
end

% try
tbl2 = readtable('tmp_output_A.txt','Delimiter','\t');
listfields = fieldnames(tbl2);
if(isScaled==3)
  for ll=2:(length(listfields)-1)
    tmp_centrality =(getfield(tbl2,listfields{ll}));
    %tmp_centrality = zscore(tmp_centrality);
    tmp_centrality = tmp_centrality/sum(tmp_centrality);
    tbl2 = setfield(tbl2,listfields{ll},tmp_centrality);
  end
end

centrality = [];
for ll=2:(length(listfields)-1)
	centrality = cat(2,centrality,table2array(tbl2(:,listfields{ll})));
end

% centrality = table2array(tbl2(:,'betweennessCentrality'));
% centrality = cat(2,centrality,table2array(tbl2(:,'closenessCentrality')));
% centrality = cat(2,centrality,table2array(tbl2(:,'closenessVitality')));

unix('rm tmp_A.txt tmp_output_A.txt');


end
