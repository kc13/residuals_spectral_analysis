function [stat_info] = get_spec_stats(pow,f,srch_bnds,cntl_bnds,Alpha)
%get_spec_stats.m  5/8/22
% this will return several fields
% representing one- and two-tailed tests
% and outputs that are useful different purposes (e.g., plotting)

    %confirm same orientation
    f = reshape(f,size(pow));
	cntl_inds = find(f > cntl_bnds(1) & f <= cntl_bnds(2));
	f_ix = f > srch_bnds(1) & f <= srch_bnds(2);
	srch_inds = find(f_ix);
	ff = f(srch_inds);
	AlphaC1T = (Alpha/length(srch_inds)); % 1-tailed
	AlphaC2T = AlphaC1T/2;
	sig_thresh_1T = norminv(1-AlphaC1T)*std(pow(cntl_inds))+mean(pow(cntl_inds));
	sig_thresh_2T = norminv(1-AlphaC2T)*std(pow(cntl_inds))+mean(pow(cntl_inds));
	%% note these vectors limited to search region length
	mnF = pow(f_ix);
	sig_ix_1T = mnF > sig_thresh_1T;
	sig_ix_2T = mnF > sig_thresh_2T;
    sig_ix_1TL = (pow > sig_thresh_1T) & f_ix;
    sig_ix_2TL = (pow > sig_thresh_2T) & f_ix;
    sigvec1T = nan(size(mnF));
    sigvec2T = nan(size(mnF));
    sigvec1T(sig_ix_1T) = mnF(sig_ix_1T);
    sigvec2T(sig_ix_2T) = mnF(sig_ix_2T);
    sigvec1TL = nan(size(pow));
    sigvec2TL = nan(size(pow));
    sigvec1TL(sig_ix_1TL) = pow(sig_ix_1TL);
    sigvec2TL(sig_ix_2TL) = pow(sig_ix_2TL);
    
    % pack into struct
    stat_info.f_ix = f_ix;
    stat_info.ff = ff;
    stat_info.mnF = mnF;
    stat_info.pow = pow;
    stat_info.f = f;
    stat_info.AlphaC1T = AlphaC1T; 
    stat_info.AlphaC2T = AlphaC2T; 
    stat_info.sig_thresh_1T = sig_thresh_1T;
    stat_info.sig_thresh_2T = sig_thresh_2T;
    stat_info.sig_ix_1T = sig_ix_1T;
    stat_info.sig_ix_2T = sig_ix_2T;
    stat_info.sig_ix_1TL = sig_ix_1TL;
    stat_info.sig_ix_2TL = sig_ix_2TL;
    stat_info.sigvec1T = sigvec1T;
    stat_info.sigvec2T = sigvec2T;
    stat_info.sigvec1TL = sigvec1TL;
    stat_info.sigvec2TL = sigvec2TL;    
    
end %fn