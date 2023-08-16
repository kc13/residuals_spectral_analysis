function [ax] = plotSpikeSample(data,msecs,category)
%plotSpikeSample.m 8/31/22

tickdir = 'out';
dobox = 'off';
dogrid = false;

xUB = max(msecs);
xUBR = round(xUB,-2);
xt = [0:50:xUBR]; 

clrsMap = containers.Map();
clrsMap('orig') = [1 1 1]*0.6; 
clrsMap('fit') = [0.3961 0.1098 0.0863];
clrsMap('res') = [0.0039 0.1255 0.2353]; 
titleMap = containers.Map();
titleMap('orig') = 'data segment';
titleMap('fit') = 'GLM fit'; 
titleMap('res') = 'GLM residuals';

xlblMap = containers.Map();
xlblMap('orig') = '';
xlblMap('fit') = 'timestamp (1 ms bins)';
xlblMap('res') = '';
ylblMap = containers.Map();
ylblMap('orig') = 'zoomed (max = 1)';
ylblMap('fit') = '';
ylblMap('res') = '';

s = stem(msecs,data(msecs),'Color',clrsMap(category));
s.Marker = 'o';
s.ShowBaseLine = 'off';
s.MarkerFaceColor = s.Color;
s.MarkerEdgeColor = 'none';
s.MarkerSize = 1;
s.LineWidth = 0.5;

ax = gca;
ax.TickDir = tickdir;
ax.Box = dobox;
ax.XLim = [0 max(msecs)];  
xlabel(xlblMap(category))
ylabel(ylblMap(category))
title(titleMap(category))
if dogrid; grid on; end
   
end