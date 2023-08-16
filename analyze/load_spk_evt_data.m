function [S_out] = load_spk_evt_data(Recordings, datadir, inds, rt_good, mt_good)
%load_spk_evt_data.m
% reflects code contributions from RST, DK, KC

t_fields = {'starttime' 'endtime' 'rewardon' 'mvt_onset' ...
	'mvt_end' 'return_onset' 'cue_onset' 'return_end' 'rewardoff'};

chan = Recordings.Chan(inds(1));
sort = Recordings.Sort(inds(1));
area = Recordings.Area(inds(1));
sort_qual = Recordings.Qualitypostsorting(inds(1));
if iscell(sort_qual)
	sort_qual = sort_qual{1};
end

file_start_time = 0;
CumEvts = ([]);

for i=1:length(inds)
	ind = inds(i);
	fname = Recordings.Filename{ind};
	uname = [fname '_' num2str(chan) '_' num2str(sort)];
	disp(['Loading ', fname, ' Session:',num2str(Recordings.Session(ind)),' Ch:',num2str(Recordings.Chan(ind)),' Sort:',num2str(Recordings.Sort(ind)),' RecInds= ', num2str(inds(1))]);
	
	filepath = fullfile(datadir,fname);
	load(filepath);
	
	Evts = S.Evts;
	session_end = max(Evts.endtime);
	
	RT = Evts.mvt_onset - Evts.cue_onset;
	MT = Evts.mvt_end - Evts.mvt_onset;
	
	% good period applies if acceptable isolation 
	% of the unit was only achieved during one or more
	% specific intervals of the recording
	if ~isempty(Recordings.GoodPeriod{ind})
		strs = strsplit(Recordings.GoodPeriod{ind},',');
		nSt = numel(strs);
		gp = cell(nSt,1);
		for st = 1:nSt
			str = strs{st};
			if strfind(str,'end')
				goodper(1) = sscanf(str,'%i - end')';
				goodper(2) = max(Evts.endtime);
			else
				goodper = sscanf(str,'%i - %i')';
			end %if
			gp{st} = goodper;
		end %for st
		gpt = cell2mat(gp)'; 
		good_trials = find(~isnan(Evts.rewardon) & inrange(RT,rt_good) & inrange(MT,mt_good) & ...
			any(Evts.starttime>(gpt(1,:)) & Evts.endtime<(gpt(2,:)),2))';		
	else 
		good_trials = find(~isnan(Evts.rewardon) & inrange(RT,rt_good) & inrange(MT,mt_good))';		
		gp = {[0,max(Evts.endtime)]}; 
	end %if ~isempty gp
	
	% From each field of Evts structure
	Evtnames = fieldnames(Evts);
	for f=1:length(Evtnames)
		% Select good trials 
		tmp = Evts.(Evtnames{f})(good_trials) ;
		% Add file_start_time to any time-coding fields of Evts
		if any(strcmp( Evtnames{f}, t_fields))
			if ~isfield(CumEvts,Evtnames{f})
				CumEvts.(Evtnames{f}) = (tmp + file_start_time);
			else
				CumEvts.(Evtnames{f}) = [CumEvts.(Evtnames{f}) ; (tmp + file_start_time) ];
			end
		else  % Don't add file_start_time to non-time fields
			if ~isfield(CumEvts,Evtnames{f})
				CumEvts.(Evtnames{f}) = tmp;
			else
				CumEvts.(Evtnames{f}) = [CumEvts.(Evtnames{f}) ; tmp ];
			end
		end
	end	
	
	% field to code which session a trial is from
	tmp = repmat(i,length(good_trials),1);
	if ~isfield(CumEvts,'Session')
		CumEvts.Session = tmp ;
	else
		CumEvts.Session = [CumEvts.Session ; tmp];
	end	
	
	%trials within session: not incremented cumulatively
	tmp = [1:length(Evts.trialnum)]';
	trialcode = tmp(good_trials);
	if ~isfield(CumEvts,'trialcode')
		CumEvts.trialcode = trialcode;
	else
		CumEvts.trialcode = [CumEvts.trialcode; trialcode];
	end

	% filename field
	tmp = repmat(fname,length(good_trials),1);
	if ~isfield(CumEvts,'Filename')
		CumEvts.Filename = tmp ;
	else
		CumEvts.Filename = [CumEvts.Filename ; tmp];
	end   	

	% trial IDs
	trstr = cell2mat(arrayfun(@(x) sprintf('%03u',x),trialcode,'UniformOutput',false));	
	trID = cellstr(horzcat(tmp,repmat('_TC_',numel(good_trials),1),trstr)); 
	if ~isfield(CumEvts,'trID')
		CumEvts.trID = trID ;
	else
		CumEvts.trID = [CumEvts.trID ; trID];
	end  	


	clear kch kst
	knu = length(S.units);
	for kn=1:knu
		kch(kn) = S.units(kn).chan;	% Chan #
        if ischar(S.units(kn).sort)
    		kst(kn) = (double(S.units(kn).sort)-96);	% Sort #
        else
            kst(kn) = S.units(kn).sort;
        end
	end

	Uind = find( kch == chan & kst == sort );
	if isempty(Uind)
		error('Cannot find correct unit #!')
	elseif length(Uind)>1
		error('Found too many matching units!')
    end

    S.units(Uind).uname=[fname '_' num2str(S.units(Uind).chan) '_' num2str(S.units(Uind).sort)];

	if ~exist('units','var')
		units = S.units(Uind);
		units.brain_area = area;
		units.fname = fname;
        units.uname = uname;
		units.sort_qual = sort_qual;
        units.gp = gp;
        units.Uind = Uind;
	else
		units.fname = [units.fname '/' fname];
		units.uname = [units.uname '/' S.units(Uind).uname];
		units.ts = [units.ts; (S.units(Uind).ts + file_start_time) ];
        temp_gp = cell(length(gp),1);
        for g = 1:length(gp)
            temp_gp{g,1} = gp{g} + file_start_time;
        end
        units.gp = [units.gp;temp_gp];
        units.Uind = [units.Uind; Uind];
	end
	
	file_start_time = file_start_time + session_end + 10;
end %for

if exist('units','var')
	S_out.units = units;
	S_out.Evts = CumEvts;
	S_out.info = S.info;  % this will have fs
else
	S_out = ([]);
end % if


end %fn

	