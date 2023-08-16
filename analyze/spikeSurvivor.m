function [S,lags] = spikeSurvivor(delta_mat,time_dim)
%spikeSurvivor.m
% inputs: 
% delta_mat: binary matrix, spike occurrence = 1
% time_dim: 1 if time down rows, 2 if time across cols
% outputs: 
% S: vector (row, col for time_dim 2,1)
% survival function up through max lag in delta_mat
% lags: survival times at which S is evaluated
% using discrete survival convention of ISI >= isi(j)
% see https://data.princeton.edu/wws509/notes/c7s6
% so S(1) = 100%


% time<->row convention
if time_dim == 1
    dmat = delta_mat';
else
    dmat = delta_mat;
end %if

ISI = getISI(dmat);
edges = [1:max(ISI)+1];
cts = histcounts(ISI,edges);
S = cumsum(cts,'reverse')./sum(cts);
lags = edges(1:end-1);
% >= convention
% see https://data.princeton.edu/wws509/notes/c7s6

if time_dim == 1
    S = S';
    lags = lags';
end
    
end %fn

