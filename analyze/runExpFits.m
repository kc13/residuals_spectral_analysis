function [mdl,dt,p] = runExpFits(isi_pdf,isiLags,run_parfor)
%runExpFits.m

if ~exist('run_parfor','var')
 run_parfor = true;
end   

n = numel(isiLags);
mdl = cell(n-1,1);
dt = cell(n-1,1); %outputs tables
p = nan(n-1,1);
disp('running ISI pdf fits')

if run_parfor
	poolobj = gcp; % start if not started yet
    parfor m = 1:n-1
        Y = isi_pdf(m:end)/sum(isi_pdf(m:end)); % renorm
        X = 1:length(Y);
        [mdl{m}] = fitglm(X,Y,'linear','Distribution','poisson');
        dt{m} = devianceTest(mdl{m});
        p(m) = dt{m}.pValue(2);
    end %m	
else % not parfor
    for m = 1:n-1
        Y = isi_pdf(m:end)/sum(isi_pdf(m:end)); % renorm
        X = 1:length(Y);
        [mdl{m}] = fitglm(X,Y,'linear','Distribution','poisson');
        dt{m} = devianceTest(mdl{m});
        p(m) = dt{m}.pValue(2);
    end %m 	
end %if

end %fn
