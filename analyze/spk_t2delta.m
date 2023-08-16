function spk_delta = spk_t2delta(spk_t,data_len,msec)
% Function to produce a delta function from a list of spike times
%
%	function spk_delta = spk_t2delta(spk_t,data_len,msec)
%
% Enter as arguement:
%
%       t_spk = times of spikes in seconds - may be a vector or array
%			arrays assume rows = trials
%       data_len = length of delta function to be returned (in msec - optional)
%		msec = flag TRUE indicates that input is in msec (optional, default
%			is seconds)
%	msec = optional flag to indicate spike times are in milliseconds
%
% Returns:
%       spk_delta = msec resolution delta function 
%                   (spike = 1, no-spike = 0)
%
%   RST March 25, 2002
%   	Oct 03, allow option of no data len
%	RST Nov, 2003	allow array of spike times
%	RST Aug, 2005 check for and correct negative spike times
%
seconds = true;

% Remove negative offset if present
mn = min(min(spk_t));
if mn<0
	spk_t = spk_t -mn;
end

% Discover data length (in msec) if not specified
if ~exist('data_len','var')
	data_len = 1+round(1000*max(nanmax(spk_t,[],2)));	
end

if exist('msec','var')
	if msec,	seconds = false;	end
end
% if seconds
% 	data_len = round(data_len*1000);
% end

[n_row,width] = size(spk_t);
% Better to not try to predict if the data is correctly oriented
% rows = trials
% Edited by HW 181024
% if width < 2 && n_row > 1	% Rotate if times in column
%  	spk_t = spk_t';
%  	[n_row,width] = size(spk_t);
% end
	
spk_delta = zeros(n_row,data_len); % make delta fn as long as required

for i = 1:n_row
	% Convert spike train to delta function
	if seconds
		tmp_spk = 1+round(1000*spk_t(i,:));        % put spk time in msec starting @ 1 msec
	else
		tmp_spk = 1+spk_t(i,:);
	end
	tmp_spk(isnan(tmp_spk)) = [];		% delete NaNs
	spk_delta(i,squeeze(tmp_spk)) = 1;		% code times of spikes as 1's
end
	

