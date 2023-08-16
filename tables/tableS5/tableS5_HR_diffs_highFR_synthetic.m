%tableS5_HR_diffs_highFR_synthetic.m
%%
clear all;  
close all;
clc;
warning('on','all')

scriptdirs = {pwd; '..\..\plot'; '..\..\analyze'; '..\..\helper_functions'};
addpath(scriptdirs{:})

datadir = pwd;
ssfile = 'subsample1000_highFR_stdRP.mat';   
sspath = fullfile(datadir,ssfile);

%%

tabledir = datadir;
tablefile = 'HRLM_3I_lag_highFR.xlsx'; 
tablepath = fullfile(tabledir,tablefile);
%%
writeTbl = true;

%%
fprintf('loading %s\n',sspath)
load(sspath,'T');

%%
T.HRDRC = -1*T.HRDCR;
T.FADRC = -1*T.FADCR;

%% exclude m = 0
% because HR is zero by definition
ix = T.M ~= 0;
yh = T.HRDRC(ix);
nw = T.NW(ix);
fo = T.FO(ix);
pbo = T.PBO(ix);
m = T.M(ix); 
%%
dmn = @(x) x-mean(x); 
dnw = dmn(nw);
dfo = dmn(fo);
dpbo = dmn(pbo);
dm = dmn(m);

%%
X = [dnw dfo dpbo dm dm.^2];  
tbl = array2table([X,yh]);
tbl.Properties.VariableNames = {'T','f_osc','pbOff','m','m2','HRdiff'};
vars = tbl.Properties.VariableNames;
xvars = vars(1:end-1);
exprME = ['1+ ',char(join(xvars,' + '))];
expr2I = char(join(join(nchoosek(xvars,2),':',2)',' + '));
expr3I = char(join(join(nchoosek(xvars,3),':',2)',' + '));
expr = ['HRdiff ~ ',exprME,' + ',expr2I,' + ',expr3I];

%%
mdl3I = fitlm(tbl,expr); 

if writeTbl
    toG = @(x) arrayfun(@(y) sprintf('%0.2g',y),x,'UniformOutput',false);
    toF = @(x) arrayfun(@(y) sprintf('%0.2f',y),x,'UniformOutput',false);
    toD = @(x) regexprep(x,'^0$','--');
    old = {'\(Intercept\)',':','pbOff'};
    new = {'Intercept',' x ','pbase_offset'};
    lblRep = @(x) regexprep(x,old,new);

    mdlTbl = mdl3I.Coefficients;
    newMat = [lblRep(mdlTbl.Properties.RowNames),toG(mdlTbl.Estimate),toG(mdlTbl.SE),toF(mdlTbl.tStat),toD(toG(mdlTbl.pValue))];
    vnames = {'Effect','Beta','SE','t','p'};
    newTbl = array2table(newMat,'VariableNames',vnames);
    
    fprintf('saving %s\n',tablepath)
    writetable(newTbl,tablepath)
end
%%