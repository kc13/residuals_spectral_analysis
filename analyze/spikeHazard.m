function [H,lags] = spikeHazard(delta_mat,time_dim)
%spikeHazard.m 5/20/22
% inputs: 
% delta_mat: binary matrix, spike occurrence = 1
% time_dim: 1 if time down rows, 2 if time across cols
% outputs: 
% H: vector (row, col for time_dim 2,1)
% hazard function up through max lag in delta_mat
% lags: using discrete survival convention of ISI >= isi(j)
% see https://data.princeton.edu/wws509/notes/c7s6
% and spikeSurvival.m


% time<->row convention
if time_dim == 1
    dmat = delta_mat';
else
    dmat = delta_mat;
end %if

ISI = getISI(dmat);
edges = [1:max(ISI)+1];
cts = histcounts(ISI,edges);
PDF = cts./sum(cts);
[S,lags] = spikeSurvivor(delta_mat,time_dim);
H = PDF./S;


end

