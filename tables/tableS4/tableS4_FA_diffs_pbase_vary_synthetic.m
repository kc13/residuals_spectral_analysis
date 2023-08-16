%tableS4_FA_diffs_pbase_vary_synthetic.m
%%
clear all;  
close all;
clc;
warning('on','all')

scriptdirs = {pwd; '..\..\plot'; '..\..\analyze'; '..\..\helper_functions'};
addpath(scriptdirs{:})

datadir = pwd;
ssfile = 'subsample1000_pbase_vary.mat';
sspath = fullfile(datadir,ssfile);


%% paths
tabledir = datadir;
tablefile = 'FALM_3I_lag_follow_up.xlsx';
tablepath = fullfile(tabledir,tablefile);
%%

%%
writeTbl = true; 

%%
fprintf('loading %s\n',sspath)
load(sspath,'T')


%% flip col
T.HRDRC = -1*T.HRDCR;
T.FADRC = -1*T.FADCR;

%% exclude m = 0
% because HR is zero by definition
% and applying same to FA to be consistent
ix = T.M ~= 0;
yh = T.FADRC(ix);
nw = T.NW(ix);
pb = T.PB(ix);
pbo = T.PBO(ix);
m = T.M(ix); 
%%
dmn = @(x) x-mean(x); 
dnw = dmn(nw);
dpb = dmn(pb);
dpbo = dmn(pbo);
dm = dmn(m);

%%
X = [dnw dpb dpbo dm dm.^2];  
tbl = array2table([X,yh]);
tbl.Properties.VariableNames = {'T','pbase','pbOff','m','m2','FAdiff'};
vars = tbl.Properties.VariableNames;
xvars = vars(1:end-1);
exprME = ['1+ ',char(join(xvars,' + '))];
expr2I = char(join(join(nchoosek(xvars,2),':',2)',' + '));
expr3I = char(join(join(nchoosek(xvars,3),':',2)',' + '));
expr = ['FAdiff ~ ',exprME,' + ',expr2I,' + ',expr3I];

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