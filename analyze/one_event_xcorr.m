function [xc_out] = one_event_xcorr(src_spk,targ_spk,stevt,endevt,xcorr_prefs)
% one_event_xcorr.m
% for autocorr, src = targ

assert(isvector(src_spk) && isvector(targ_spk),'only vector inputs accepted for spks');

store_input = optCheck(xcorr_prefs,'store_input',false,[true false]);
lags_msec = optCheck(xcorr_prefs,'lags_msec',200,[]);

% to ensure col vec (can assume input is vector)
src_spk_c = reshape(src_spk,[],1);
targ_spk_c = reshape(targ_spk,[],1);

len_msec = (endevt - stevt)*1000;

% note currently inclusive lb and exclusive ub
src_raster = one_event_raster(src_spk,stevt,endevt);
targ_raster = one_event_raster(targ_spk,stevt,endevt);

src_delta = spk_t2delta_vector(src_raster,len_msec);
targ_delta = spk_t2delta_vector(targ_raster,len_msec);

%pre-allocate
xc = nan((2*lags_msec)+1,1);
be = xc;
xcraw = xc;
xcn = xc;
xc_out.n_trials = 1;

targ_FR = sum([targ_delta(:)]/(numel([targ_delta(:)])/1000));
src_FR = sum([src_delta(:)]/(numel([src_delta(:)])/1000));
xc_out.targ_FR = targ_FR;
xc_out.src_FR = src_FR;

if store_input
    xc_out.src_delta = src_delta;
    xc_out.targ_delta = targ_delta;
end

[xc,be,xcraw] = getCCF(src_delta,targ_delta,lags_msec);

xc_out.xc = xc;
xc_out.be = be;
xc_out.xcraw = xcraw;

xcSraw = xcorr(src_delta,src_delta,lags_msec,'none');
xcTraw = xcorr(targ_delta,targ_delta,lags_msec,'none');

% like matlab scaleopt = 'normalized'
ix0 = lags_msec+1;
normterm = sqrt(xcSraw(ix0)*xcTraw(ix0));
xc_out.xcn = xcraw./normterm;
xc_out.xcn2 = despkACF(xc_out.xcn);

xc_out.nTargSpk = sum([targ_delta(:)]);
xc_out.nSrcSpk = sum([src_delta(:)]);
xc_out.denom = (targ_FR*sum([src_delta(:)]));

end %fn
