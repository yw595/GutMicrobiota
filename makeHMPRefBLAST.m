extraDir = '/mnt/extra/blast/HMPRef';
fileNames = dir(extraDir);
generaToFileNames = containers.Map;
for i=1:length(fileNames)
    % INSIDIOUS BUG, without underscore, will match genusName.faa
    % from previous runs, then genusName variable will be empty,
    % resulting in bunch of fileNames linked to ''
    if ~isempty(regexp(fileNames(i).name,'_.*\.faa$'))
        underscoreIdx = regexp(fileNames(i).name,'_');
        genusName = fileNames(i).name(1:underscoreIdx-1);
        %disp(genusName)
        if isKey(generaToFileNames,genusName)
            genusFileNames = generaToFileNames(genusName);
            genusFileNames{end+1} = fileNames(i).name;
            generaToFileNames(genusName) = genusFileNames;
        else
            generaToFileNames(genusName) = {fileNames(i).name};
        end
    end
end

if 0
generaNames = keys(generaToFileNames);
for j=1:length(generaNames)
    genusName = generaNames{j};
    disp(genusName)
    catCmd = 'cat ';
    genusFileNames = generaToFileNames(genusName);
    for i=1:length(genusFileNames)
        catCmd = [catCmd extraDir filesep genusFileNames{i} ' '];
    end
    catCmd = [catCmd '> ' extraDir filesep 'temp.faa'];
    status=system(catCmd);
    sedCmd = sprintf('%s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/modHeader.sh'],[extraDir filesep 'temp.faa'],[extraDir filesep genusName '.faa']);
    status=system(sedCmd);
    status=system(['rm ' extraDir filesep 'temp.faa']);

    remDupCmd = sprintf('%s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/removeDupSeqs.sh'],[extraDir filesep genusName '.faa'],[extraDir filesep genusName '2.faa']);
    status=system(remDupCmd);

    trimCmd = sprintf('%s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/trimLines.sh'],[extraDir filesep genusName '2.faa'],[extraDir filesep 'temp.faa']);
    status=system(trimCmd);
    
    dbCmd = ['makeblastdb -in ' extraDir filesep genusName '2.faa -dbtype prot -out ' extraDir filesep genusName ' -parse_seqids -title ' genusName];
    status=system(dbCmd);

    %status = system(['rm ' extraDir filesep genusName '.faa']);
end
end

tempFileNames = dir(tempDir);
for i=1:length(tempFileNames)
    if ~isempty(regexp(tempFileNames(i).name,'.fasta'))
        genusName = tempFileNames(i).name(1:regexp(tempFileNames(i).name,'.test.fasta')-1);
        if ~strcmp(genusName,'Shigella')
            sedCmd = sprintf('%s %s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/modHeader.sh'],[tempDir filesep genusName '.test.fasta'],[tempDir filesep genusName '.faa'],'true');
            status=system(sedCmd);
        else
            sedCmd = sprintf('%s %s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/modHeader.sh'],[tempDir filesep genusName '.test.fasta'],[tempDir filesep genusName '.faa'],'false');
            status=system(sedCmd);
            status=system(sprintf('cut -d"[" -f1 < %s > %s',[tempDir filesep genusName '.faa'],[tempDir filesep 'temp.faa']));
            status=system(sprintf('mv %s %s',[tempDir filesep 'temp.faa'],[tempDir filesep genusName '.faa']));
        end
        %status=system(sedCmd);

        remDupCmd = sprintf('%s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/removeDupSeqs.sh'],[tempDir filesep genusName '.faa'],[tempDir filesep genusName '2.faa']);
        status=system(remDupCmd);

        trimCmd = sprintf('%s %s %s',['/home/ubuntu/MATLAB/GutMicrobiota/trimLines.sh'],[tempDir filesep genusName '2.faa'],[tempDir filesep 'temp.faa']);
        status=system(trimCmd);
        
        dbCmd = ['makeblastdb -in ' tempDir filesep genusName '2.faa -dbtype prot -out ' tempDir filesep genusName ' -parse_seqids -title ' genusName];
        status=system(dbCmd);
    end
end