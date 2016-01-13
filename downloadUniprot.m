outputDir1 = '/home/ubuntu/MATLAB/GutMicrobiota/output/downloadUniprot';
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end

files = {'supernumerary.txt','butyrate.txt','folate.txt'};

for j=2:length(files)

    FI = fopen(['/home/ubuntu/MATLAB/GutMicrobiota/input' filesep files{j}]);
    dataFields = textscan(FI, repmat('%s',1,2),'Delimiter','\t', 'HeaderLines',0);
    rawTable = [dataFields{:}];
    names = rawTable(1:end,1);
    ids = rawTable(1:end,2);
    fclose(FI);
    
    for i=1:length(ids)
        status=system(sprintf('wget -O %s www.uniprot.org/uniprot/%s.fasta',[outputDir1 filesep names{i} '.faa'],ids{i}));
    end
end