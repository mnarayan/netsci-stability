% Network Metrics: Either BCT or matlab-bgl or custom

METRICS_PATH= '../packages/BCT'
addpath(genpath(METRICS_PATH));
addpath(('python'))
addpath(('netfun'))

% visualization
PLOTLY_PATH = '../packages/MATLAB-api'
addpath(genpath(PLOTLY_PATH));
addpath(genpath(fullfile(matlabroot,'toolbox','plotly')),'-end');

echo off;
warning off;
