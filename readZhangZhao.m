cd /home/ubuntu/MATLAB/GutMicrobiota
topDir = '/home/ubuntu/MATLAB/GutMicrobiota';

FI = fopen(['input' filesep 'srep14405-s2-mod.txt']);
dataFields = textscan(FI, ...
    repmat('%s',1,16),'Delimiter','\t','HeaderLines',2);
fclose(FI);
rawTable = [dataFields{:}];
% skip bottom two lines with *compared with HFD comment
rawTable = rawTable(1:end-2,:);
taxonomies = rawTable(1:end,3);
% strsplit on ; always leads to space at end if no genus, so flag on this
hasGeneraTemp = cellfun(@(x) strsplitYiping(x,';'), taxonomies, 'UniformOutput',0);
hasGenera = cellfun(@(x) ~strcmp(' ',x(5)), hasGeneraTemp);
% UniformOutput==1 gives one cell array, else cell array of 1x1 cell
% arrays???, and fail???
generaTemp = cellfun(@(x) x(end), hasGeneraTemp(hasGenera), 'UniformOutput', 1);
% remove beginning whitespace
genera = cellfun(@(x) x(2:end), generaTemp, 'UniformOutput', 0);
contVsB200 = rawTable(hasGenera,9);
contVsB100 = rawTable(hasGenera,10);
contVsM200 = rawTable(hasGenera,11);
ZhangZhaoUpGenera = unique(genera(strcmp(contVsB200,'up') & strcmp(contVsB100,'up') & strcmp(contVsM200,'up')));
ZhangZhaoDownGenera = unique(genera(strcmp(contVsB200,'down') & strcmp(contVsB100,'down') & strcmp(contVsM200,'down')));
ZhangZhaoDownGenera = changeNames(ZhangZhaoDownGenera);
ZhangZhaoUpGenera = changeNames(ZhangZhaoUpGenera);
ZhangZhaoDownStrictGenera = ZhangZhaoDownGenera(~ismember(ZhangZhaoDownGenera,ZhangZhaoUpGenera));
ZhangZhaoUpStrictGenera = ZhangZhaoUpGenera(~ismember(ZhangZhaoUpGenera,ZhangZhaoDownGenera));