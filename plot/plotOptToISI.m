function [ax] = plotOptToISI(rpInfo,priColors,figprefs)
%plotOptToISI.m 

if isempty(priColors)
    priColors = containers.Map();
    priColors('HL') = 'k';
    priColors('bars') = [0.6 0.6 0.6];
end


isiLags = rpInfo.isiLags;
isiPDF = rpInfo.isiPDF;

tickdir = 'out';
dobox = 'off';
dogrid = false;
%9/1/22
doleg = optCheck(figprefs,'doleg',true,[false true]);

b = bar(isiLags,isiPDF,'FaceColor',priColors('bars'));
hold all;
if dogrid; grid on; end

optMdl = rpInfo.optMdl;
HLclrs = priColors('HL');
HLix = rpInfo.optLag;
nHL = numel(HLix);
goArr = gobjects(nHL,1);
lw = 1; %0.5 default

for h = 1:nHL
    lag = rpInfo.optLag;
    mdl = optMdl;
    mdlX = mdl.Variables.x1; 
    YHat = predict(mdl,mdlX);  % mdl handles intercept
    YHatScale = YHat*(sum(isiPDF(lag:end)));
    XShift = mdlX + lag -1;
    goArr(h) = plot(XShift,YHatScale,'Color',HLclrs(h,:),'LineWidth',lw);
end %h

if isfield(figprefs,'xlbl')
    xlabel(figprefs.xlbl)
else
    xlabel('lags (ms)')
end
ylabel('p(ISI = lag)')

if doleg
    legtext = regexprep(cellstrnum(toCol(HLix)),'(.+)','start lag = $1 ms');  
    legend(goArr,legtext,'Box','off')
end

hold off;
if isfield(figprefs,'xlims')
   xlim(figprefs.xlims) 
end
ylims = optCheck(figprefs,'ylims',[],[]);
if ~isempty(ylims)
    ylim(figprefs.ylims);
end

ax = gca;
ax.TickDir = tickdir;
ax.Box = dobox;

    

end %fn

