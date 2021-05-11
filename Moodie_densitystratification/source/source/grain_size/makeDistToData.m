function [long, lam] = makeDistToData(gsDist)
    gsClassCt = floor( table2array(gsDist).*10 );
    [gsClassNum] = gsDistClass2num(gsDist);
    nR = nansum( gsClassCt ); % number of rows to fill out
    nC = size(gsDist, 1); % number of grain classes
    long = NaN(nR, 1);
    f = 1; % initialize filled to idx
    for j = 2:nC
        nF = gsClassCt(j); % number of slots to fill with this grain class
        long(f:(f+nF-1)) = repmat(gsClassNum(j), nF, 1);
        f = f+nF;
    end
    longPhi = -1 .* log2(long./1000);
    longStd = std(longPhi, 'omitnan');
    lam = 1 - (0.28 .* longStd);
end