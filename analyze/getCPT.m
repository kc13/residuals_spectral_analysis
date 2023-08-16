function [c_pt] = getCPT(lag,spline_method,brkpt)
% ref for splines: https://github.com/MehradSm/Modified-Spline-Regression

nSeg = ceil(lag/brkpt);
c_pt_ctr = prctile([0:lag+1],([0:nSeg]./nSeg)*100);
switch spline_method
    case 'modcard'
       c_pt = c_pt_ctr;     
    case 'card'  
	   c_pt = [c_pt_ctr(1)-range(c_pt_ctr)/nSeg c_pt_ctr c_pt_ctr(end)+range(c_pt_ctr)/nSeg];
    otherwise
        error('unrecognized spline method')
end %sw
            
end %fn

