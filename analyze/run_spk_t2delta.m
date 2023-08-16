function spk_delta = run_spk_t2delta(spk_t,data_len)
%1/30/21
%wrapper as a post-hoc fix for existing in-house spk_t2delta function
% it can return a vector that contains one extra element (due to rounding)

    runWarn = false;

	spk_delta_raw = spk_t2delta(spk_t,data_len);
    if isvector(spk_delta_raw) 
        spk_delta = spk_delta_raw(1:data_len);
    else % spk_t2delta will return trials along rows
        spk_delta = spk_delta_raw(:,1:data_len);
    end
	if length(spk_delta) ~= length(spk_delta_raw) && runWarn
		warning('it was necessary to correct the length of the spk_t2delta output');
	end

end %fn