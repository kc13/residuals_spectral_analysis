function [rpInfo] = runRPfinder_vector(delta_vec,run_parfor)
% will output both estimated RP and
% additional useful info

if ~exist('run_parfor','var')
    run_parfor = true;
end

assert(isvector(delta_vec),'use non-vector verison for matrix data')

[isiPDF,isiLags] = isiPDF_vector(delta_vec);
[expMdl,dt,p] = runExpFits(isiPDF,isiLags,run_parfor);
[npks,locs] = findpeaks(-p); 
if isempty(locs)
	disp('no dt peak -- going into debug mode')
	keyboard
end %if
[topPk,topPkIx] = max(npks);
firstPk = npks(1);
lagIx = locs(1);
optLag = isiLags(lagIx);
RPend = optLag-1;
optMdl = expMdl{lagIx}; 
[H,hLags] = spikeHazard_vector(delta_vec);
optMdlX = optMdl.Variables.x1;
optMdlY = optMdl.Variables.y; 
[optYHat] =  predict(optMdl,optMdlX);

ISI = getISI_vector(delta_vec);

rpInfo.isiPDF = isiPDF;
rpInfo.isiLags = isiLags;
rpInfo.expMdl = expMdl;
rpInfo.dt = dt;
rpInfo.p = p;
rpInfo.npks = npks;
rpInfo.locs = locs;
rpInfo.firstPk = firstPk;
rpInfo.lagIx = lagIx;
rpInfo.optLag = optLag;
rpInfo.RPend = RPend;
rpInfo.optMdl = optMdl;
rpInfo.H = H;
rpInfo.hLags = hLags;
rpInfo.optMdlX = optMdlX;
rpInfo.optMdlY = optMdlY;
rpInfo.optYHat = optYHat;


end %fn
