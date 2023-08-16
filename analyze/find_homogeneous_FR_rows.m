function [clean_data,info] = find_homogeneous_FR_rows(data)
% find_homogeneous_FR_rows.m
% data: assumes nTr x nMS structure

minSpkTr = 2;
outlierFac = 3;

FRs = mean(data,2)*1000;
threshIX = FRs >= minSpkTr;
threshRows = find(threshIX); 

gmmData = data(threshIX,:);
gmmFRs = FRs(threshIX);
nGRows = size(gmmData,1);

clean_data = [];
cleanIX = [];
gaussIX = [];
optGMM = [];

if nGRows > 1
    maxK = nGRows;
    bicVec = nan(maxK,1);
    gmmVec = cell(maxK,1);
    for k = 1:maxK
        try
            [init_centers,init_sigma,init_props] = varPart1D(gmmFRs,k);
            S.mu = init_centers'; % k x d
            S.Sigma = nan(1,1,k); % d x d x k
            S.Sigma(1,1,:) = init_sigma; 
            S.ComponentProportion = init_props; % 1 x k
            gmmVec{k} = fitgmdist(gmmFRs,k,'Start',S);
            bicVec(k) = gmmVec{k}.BIC;
        catch
            bicVec(k) = nan;
        end
        
        if k > 1
            if bicVec(k-1) < bicVec(k) || isnan(bicVec(k))
                optK = k-1;
                optGMM = gmmVec{k-1};
                break;
            end
        end %if
    end % for k

    lbls = cluster(optGMM,gmmFRs);
    modeLbl = mode(lbls);
    gmmToModeIX = lbls == modeLbl;
    optMu = squeeze(optGMM.mu);
    modeMu = optMu(modeLbl);
    optSigma = squeeze(optGMM.Sigma);
    modeSigma = optSigma(modeLbl);
    modeSD = sqrt(modeSigma);
    modeRows = threshRows(gmmToModeIX);
    gaussIX = zeros(size(threshIX));
    gaussIX(modeRows) = 1;
    
    oneGFRs = gmmFRs(gmmToModeIX);
    oneGData = gmmData(gmmToModeIX,:);
    
    modeToNoOutIX = abs(oneGFRs-modeMu) <= outlierFac*modeSD;  
    noOutFRs = oneGFRs(modeToNoOutIX);
    noOutData = oneGData(modeToNoOutIX,:);  

    clean_data = noOutData; 
    cleanRows = modeRows(modeToNoOutIX);
    cleanIX = zeros(size(threshIX));
    cleanIX(cleanRows) = 1;
end %if

info.threshIX = threshIX;
info.cleanIX = cleanIX;
info.gaussIX = gaussIX;
info.nClean = nnz(cleanIX);
info.optGMM = optGMM;
info.nStart = numel(FRs);

end