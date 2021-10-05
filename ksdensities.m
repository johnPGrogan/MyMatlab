function ksdensities(yData,varargin)
% loop through each column and do a ksdensity plot
% hold on
% passes varargin to each

%%
nCols = size(yData,2);

[fout0,xout,u,ksinfo] = deal(cell(nCols,1));
for i = 1:nCols
    hold on;
    ksdensity(yData(:,i),varargin);
end