topDir = '/home/ubuntu/MATLAB/GutMicrobiota';
FI = fopen([topDir filesep 'input' filesep 'HMRGD.csv']);
dataFields = textscan(FI, ...
    repmat('%s',1,18),'Delimiter',',','HeaderLines',1);
fclose(FI);
rawTable = [dataFields{:}];
%get rid of double quotes surround words
species = cellfun(@(x) x(2:end-1),rawTable(1:end,1),'UniformOutput',0);
locations = cellfun(@(x) x(2:end-1),rawTable(1:end,2),'UniformOutput',0);
isValid = ~strcmp(species,'') & strcmp(locations,'Gastrointestinal_tract');
rawHMPRefSpecies = species(isValid);
HMPRefSpecies = rawHMPRefSpecies;
for i=1:length(HMPRefSpecies)
    HMPRefSpecies{i} = strrep(HMPRefSpecies{i},' ','_');
    HMPRefSpecies{i} = strrep(HMPRefSpecies{i},'\.','');
    HMPRefSpecies{i} = strrep(HMPRefSpecies{i},'/','_');
end

extraDir = '/mnt/extra/blast/HMPRef';
if ~exist(extraDir,'dir')
    mkdir(extraDir);
end
inputFI = fopen(['/mnt/extra/blast/Gastrointestinal_tract.pep.fsa']);
filenamesToFIs = containers.Map;
line = fgetl(inputFI);
while line ~= -1
    if ~isempty(regexp(line,'>'))
        leftBrace = regexp(line,'\[.*\]$');
        rightBrace = regexp(line,'\]$');
        outputSpecies = line(leftBrace+1:rightBrace-1);
        % note there will one .faa for genes without any species
        outputFilename = [HMPRefSpecies{strcmp(rawHMPRefSpecies, outputSpecies)} '.faa'];
        % cannot use file existence, previous runs may already have made file
        if ~isKey(filenamesToFIs,outputFilename) %~exist([extraDir filesep outputFilename],'file')
            newFI = fopen([extraDir filesep outputFilename],'w');
            filenamesToFIs(outputFilename) = newFI;
            disp(outputSpecies);
        else
            newFI = filenamesToFIs(outputFilename);
        end
        line = line(1:leftBrace-2);
    end
    fprintf(newFI,'%s\n',line);
    line = fgetl(inputFI);
end

% note only ~425 files, although 455 entries in HMPRefSpecies, some
% species not found in Gastrointestinal_tract???
for newFI = values(filenamesToFIs)
    fclose(newFI);
end
fclose(inputFI);