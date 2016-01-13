FI = fopen('/mnt/extra/blast/assembly_result.txt');
dataFields = textscan(FI,'%s%s%s%s', 'Delimiter', '\t', 'HeaderLines',1);
fclose(FI);
FI = fopen('/mnt/extra/blast/refseqIDsToSpecies.txt');
refseqIDsToSpecies = containers.Map;
line = fgetl(FI);
while line~=-1
    words = strsplit(line,'\t');
    words1 = strsplit(words{2},' ');
    refseqIDsToSpecies(words{1}) = [words1{1} ' ' words1{2}];
    line = fgetl(FI);
end
fclose(FI);

rawTable = [dataFields{:}];
refseqIDs = rawTable(1:end,3);
tempIdx = 120;
for i=tempIdx:tempIdx%length(refseqIDs)
    % fullID = refseqIDs{i};
    % nonZeroDigits = fullID(regexp(fullID,'0[1-9]')+1:regexp(fullID,'\.')-1);
    % versionNumber = fullID(regexp(fullID,'\.')+1:end);
    % fullID = [fullID '_ASM' nonZeroDigits(1:end-1) 'v' versionNumber];
    % speciesName = strrep(refseqIDsToSpecies(refseqIDs{i}),' ','_');
    % % two Pseudomonas_syringae in refseqIDsToSpecies, only one in
    % % ncbi reference sequence
    % if strcmp(speciesName,'Pseudomonas_syringae')
    %     fullID = 'GCF_000012245.1_ASM1224v1';
    % end
    % if strcmp(speciesName,'Klebsiella_pneumoniae')
    %     fullID = refseqIDs{i};
    %     fullID = [fullID '_ASM' nonZeroDigits(1:end-1) 'v2'];
    % end
    % if strcmp(speciesName,'Vibrio_fischeri')
    %     speciesName = 'Aliivibrio_fischeri';
    % end
    % if strcmp(speciesName,'Caulobacter_crescentus')
    %     speciesName = 'Caulobacter_vibrioides';
    % end
    % if strcmp(speciesName,'Chlorobium_tepidum')
    %     speciesName = 'Chlorobaculum_tepidum';
    % end
    % if strcmp(speciesName,'Chlamydophila_pneumoniae')
    %     speciesName = 'Chlamydia_pneumoniae';
    % end
    % if strcmp(speciesName,'Lactobacillus_casei')
    %     continue;
    % end
    % status=system(['wget -O /mnt/extra/blast' filesep speciesName '.faa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria' filesep speciesName filesep 'reference' filesep fullID filesep fullID '_protein.faa.gz']);
end

dirFiles = dir('/mnt/extra/blast');
%count = 0;
for i=1:length(dirFiles)
    if ~isempty(regexp(dirFiles(i).name,'.faa$')) %&& count < 2
        %count = count+1;
        speciesName = dirFiles(i).name(1:regexp(dirFiles(i).name, '.faa')-1);
        
        disp(speciesName)
        FI = fopen(['/mnt/extra/blast' filesep speciesName '.faa']);
        dirName = ['/mnt/extra/blast' filesep speciesName 'RefSeq'];
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        line = fgetl(FI);
        isETC=0;
        count=0;
        while line~=-1
            count = count+1;
            %disp(count)
            if ~isempty(regexp(line,['(ATP synthase|NADH:ubiquin|NADH dehydrogenase|cytochrome|succinate dehydrogenase)']))
                %disp(line)
                tempLine = line(regexp(line,'\d ')+2:regexp(line,' [')- 1);
                tempLine = strrep(tempLine,'/','');
                tempLine = strrep(tempLine,' ','_');
                %disp(tempLine)
                FIout = fopen([dirName filesep tempLine '.faa'],'w');
                isETC=1;
            elseif ~isempty(regexp(line,'>')) && isETC
                isETC=0;
                fclose(FIout);
            end
            if isETC
                fprintf(FIout,'%s\n',line);
            end
            line = fgetl(FI);
        end
        fclose(FI);
    end
end