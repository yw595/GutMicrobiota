configGutMicrobiota;
runBLAST = 1;
runCDD = 0;
if runBLAST
    outputDir1 = [outputDir filesep 'runBLAST'];
elseif runCDD
    outputDir1 = [outputDir filesep 'runCDD'];
end
useNrGenomes = 0;
useHMPGenomes = 1;
count=0;
if useHMPGenomes
    outputDir1 = [outputDir1 'HMP'];
elseif useNrGenomes
    outputDir1 = [outputDir1 'Nr'];
end
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

if runBLAST
    % STUPID, but unavoidable bug only with parfor and runBLAST,
    % sees queriesAlt{z} in runCDD, thinks the index exceeds length 
    geneNamesTemp = values(complexesToGenes); geneNames = [geneNamesTemp{:}]; queries = geneNames; queriesAlt = queries;
else    
    queriesAlt = cddIDs;
    queries = cddNames;
end

%bigScores = zeros(length(complexes),1,2,2,1);
parfor z=1:length(queries)
    outputDir2 = [outputDir1 filesep queries{z}];
    if ~exist(outputDir2,'dir')
        mkdir(outputDir2)
    end
    for j=1:2
        if j==1
            upGenera = ZhangZhaoUpGenera;
            downGenera = ZhangZhaoDownGenera;
            upStrictGenera = ZhangZhaoUpStrictGenera;
            downStrictGenera = ZhangZhaoDownStrictGenera;
            outputDir3 = [outputDir2 filesep 'ZhangZhao'];
        else
            upGenera = ForslundHildebrandUpGenera;
            downGenera = ForslundHildebrandDownGenera;
            upStrictGenera = ForslundHildebrandUpStrictGenera;
            downStrictGenera = ForslundHildebrandDownStrictGenera;
            outputDir3 = [outputDir2 filesep 'ForslundHildebrand'];
        end        
        if ~exist(outputDir3,'dir')
            mkdir(outputDir3)
        end
        
        scores = [];
        for k=1:4
            if k==1
                genera = upGenera;
                outputDir4 = [outputDir3 filesep 'up'];
            elseif k==2
                genera = downGenera;
                outputDir4 = [outputDir3 filesep 'down'];
            elseif k==3
                genera = upStrictGenera;
                outputDir4 = [outputDir3 filesep 'upStrict'];
            else
                genera = downStrictGenera;
                outputDir4 = [outputDir3 filesep 'downStrict'];
            end            
            if ~exist(outputDir4,'dir')
                mkdir(outputDir4)
            end

            if k==3 || k==4
                outputDirBefore = outputDir4(1:end-6);
                for i=1:length(genera)
                    status=system(sprintf('cp %s %s',[outputDirBefore filesep genera{i} '.blast'],outputDir4));
                end
            else
                for i=1:length(genera)
                    % have to define database and genomeFile as empty,
                    % else Warning in parallel about temporary variables 
                    if runBLAST
                        database = '';
                        if useNrGenomes
                            database = 'nr';
                        elseif useHMPGenomes
                            database = [HMPDir filesep genera{i}];
                            if ~exist([database '.phr'],'file')
                                database = [tempDir filesep genera{i}];
                            end
                        end
                        query = [outputDir filesep 'downloadUniprot' filesep queries{z} '.faa'];
                        if useNrGenomes
                            gilist = [outputDir filesep 'extractGIs' filesep genera{i} '.txt'];
                        else
                            gilist = '';
                        end
                        out = [outputDir4 filesep genera{i} '.blast'];
                        runBLASTOnce(database,query,gilist,out);
                    elseif runCDD
                        genomeFile = '';
                        if useNrGenomes
                            genomeFile = [tempDir filesep genera{l} '.test.fasta'];
                        elseif useHMPGenomes
                            genomeFile = [HMPDir filesep genera{l} '.faa'];
                            if ~exist(genomeFile,'file')
                                genomeFile = [tempDir filesep genera{l} '.test.fasta'];
                            end
                        end
                        database = [CDDDir filesep 'smp' filesep queriesAlt{z}];
                        temp = [outputDir4 filesep 'temp.blast'];
                        out = [outputDir4 filesep genera{l} '.blast'];
                        runCDDOnce(database,genomeFile,temp,out);                
                    end
                end
            end
    
            for i=1:length(genera)
                FI = fopen([outputDir4 filesep genera{i} '.blast']);
                if FI ~= -1
                    line = fgetl(FI);
                    if line ~= -1
                        words = strsplit(line,'\t');
                        scores(end+1) = str2num(words{1});
                    else
                        scores(end+1) = -1;
                    end
                    fclose(FI);
                else
                    scores(end+1) = -2;
                end
            end
        end
    end
end