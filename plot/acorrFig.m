function [ax] = acorrFig(xcData,figopt)
% acorrFig.m 9/28/22

xc = xcData.xcn2;
be = xcData.be;
clr = 'k';

inputax = optCheck(figopt,'ax',[],[]);
tistr = optCheck(figopt,'tistr','ACF','');

if isempty(inputax)
	plot(be,xc,'Color',clr)
else
	plot(inputax,be,xc,'Color',clr)
end

ax = gca;
ax.TickDir = 'out'; ax.Box = 'off';
title(tistr)
ylabel('n(spk|lag)/n(spk)')   
xlabel('lags (ms)')

end %fn