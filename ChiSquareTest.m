function [p,chiSq,df] = ChiSquareTest(observedData)
%[p,chiSq,df] = ChiSquareTest(observedData)
%does chi squared test on observedData, 
%returns: p value, chi-squared value and degrees of freedom

    [r,c] = size(observedData);
    expectedData = zeros(r,c);
    for i = 1:c
        for j = 1:r
            expectedData(j,i) = nansum(observedData(:,i))/nansum(nansum(observedData))*sum(observedData(j,:));
        end
    end
    chiSq = nansum(nansum(((observedData-expectedData).^2)./expectedData));
    df = (r-1) * (c-1);
    p = 1 - chi2cdf(chiSq,df);
end