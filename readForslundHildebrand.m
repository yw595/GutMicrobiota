cd /home/ubuntu/MATLAB/GutMicrobiota
topDir = '/home/ubuntu/MATLAB/GutMicrobiota';

FI = fopen(['input' filesep 'SItableS9.txt']);
dataFields = textscan(FI, ...
    repmat('%s',1,15),'Delimiter','\t','HeaderLines',1);
fclose(FI);
rawTable = [dataFields{:}];
mPlusVsMinusPVals = rawTable(1:end,7);
filterIdxs = cellfun(@(x) ~isempty(str2num(x)), mPlusVsMinusPVals);
mPlusVsMinusPVals = rawTable(filterIdxs,7);
genera = rawTable(filterIdxs,14);
directions = rawTable(filterIdxs,5);
MetFIdxs = cellfun(@(x) x(1), strfind(directions,'Metf'));
DT2Idxs = cellfun(@(x) x(1), strfind(directions,'DT2'));
ForslundHildebrandUpGenera = unique(genera(~strcmp(genera,'?') & cellfun(@(x) str2num(x), mPlusVsMinusPVals) < .1 & (MetFIdxs < DT2Idxs)));
ForslundHildebrandDownGenera = unique(genera(~strcmp(genera,'?') & cellfun(@(x) str2num(x), mPlusVsMinusPVals) < .1 & (DT2Idxs < MetFIdxs)));
ForslundHildebrandDownGenera = changeNames(ForslundHildebrandDownGenera);
ForslundHildebrandUpGenera = changeNames(ForslundHildebrandUpGenera);
ForslundHildebrandDownStrictGenera = ForslundHildebrandDownGenera(~ismember(ForslundHildebrandDownGenera,ForslundHildebrandUpGenera));
ForslundHildebrandUpStrictGenera = ForslundHildebrandUpGenera(~ismember(ForslundHildebrandUpGenera,ForslundHildebrandDownGenera));