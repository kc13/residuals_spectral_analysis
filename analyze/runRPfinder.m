function [rpInfo] = runRPfinder(delta_mat,run_parfor)
%runRPfinder.m  5/30/22
% will output both estimated RP and
% info useful for visualizations

if ~exist('run_parfor','var')
    run_parfor = true;
end

assert(~isvector(delta_mat),'use vector verison for vector data')

[isi_pdf,isiLags] = isiPDF(delta_mat);
[expMdl,dt,p] = runExpFits(isi_pdf,isiLags,run_parfor);
[npks,locs] = findpeaks(-p); 
if isempty(locs)
	disp('no dt peak -- going into debug mode')
	keyboard
end %if
firstPk = npks(1);  
lagIx = locs(1); 
optLag = isiLags(lagIx); 
RPend = optLag-1;
optMdl = expMdl{lagIx}; 
[H,hLags] = spikeHazard(delta_mat,2);
optMdlX = optMdl.Variables.x1;
optMdlY = optMdl.Variables.y;   
[optYHat] =  predict(optMdl,optMdlX); %these will be one anchored -- will 

rpInfo.isiPDF = isi_pdf;
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
