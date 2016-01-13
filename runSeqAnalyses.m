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

otherGenera = keys(generaToFileNames);
otherGenera = otherGenera(~ismember(otherGenera,uniqGenera));

bigScores = -3*ones(length(complexes),50,2,4,20);
bigQCovS = -3*ones(length(complexes),50,2,4,20);
bigPIdent = -3*ones(length(complexes),50,2,4,20);
firstOccurrenceMatrix = zeros(length(complexes),50,2,4,20);
blastOccurrenceMatrix = zeros(length(complexes),50,length(uniqGenera));
bigMatchNames = cell(length(complexes),50,2,4,20);
parfor z1 = 1:length(complexes)
    subBigScores = -3*ones(50,2,4,20);
    subBigQCovS = -3*ones(50,2,4,20);
    subBigPIdent = -3*ones(50,2,4,20);
    subFirstOccurrenceMatrix = zeros(50,2,4,20);
    subBlastOccurrenceMatrix = zeros(50,length(uniqGenera));
    subBigMatchNames = cell(50,2,4,20);
    if runBLAST
        % STUPID, but unavoidable bug only with parfor and runBLAST,
        % sees queriesAlt{z} in runCDD, thinks the index exceeds length 
        geneNames = complexesToGenes(complexes{z1}); queries = geneNames; queriesAlt = queries;
    else    
        queriesAlt = cellfun(@(x) x(regexp(x,'_')+1:end),complexesToCDDs(complexes{z1}));
        queries = cellfun(@(x) x(1:regexp(x,'_')-1),complexesToCDDs(complexes{z1}));
    end
    
    for z=1:length(queries)
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

            for k=1:5
                if k==1
                    genera = upGenera;
                    outputDir4 = [outputDir3 filesep 'up'];
                elseif k==2
                    genera = downGenera;
                    outputDir4 = [outputDir3 filesep 'down'];
                elseif k==3
                    genera = upStrictGenera;
                    outputDir4 = [outputDir3 filesep 'upStrict'];
                elseif k==4
                    genera = downStrictGenera;
                    outputDir4 = [outputDir3 filesep 'downStrict'];
                else
                    genera = otherGenera;
                    outputDir4 = [outputDir2 filesep 'other'];
                end            
                if ~exist(outputDir4,'dir')
                    mkdir(outputDir4)
                end

                if k==3 || k==4
                    outputDirBefore = outputDir4(1:end-6);
                    for i=1:length(genera)
                        status=system(sprintf('cp %s %s',[outputDirBefore filesep genera{i} '.blast'],outputDir4));
                    end
                elseif k==1 || k==2 || (k==5 && j==1)
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

                if k~=5
                    for i=1:length(genera)
                        BLASTFile = [outputDir4 filesep genera{i} '.blast'];
                        [score matchName qCovS pIdent blastOcc] = readBLASTFile(BLASTFile);
                        % bigScores=indexPar(bigScores,score,[z1,z,j,k,i]);
                        % bigQCovS=indexPar(bigQCovS,qCovS,[z1,z,j,k,i]);
                        % bigPIdent=indexPar(bigPIdent,pIdent,[z1,z,j,k,i]);
                        % %bigMatchNames{z1,z,j,k,i}=matchName;
                        % if firstOccurrenceMatrix(z1,z,j,k,i) && blastOcc ~= blastOccurrenceMatrix(z1,z,strcmp(uniqGenera,genera{i}))
                        %     disp(sprintf('WARNING DIFFERENT BLAST OCC: %d %d %d %d %d',z1,z,j,k,i));
                        % end
                        % blastOccurrenceMatrix(z1,z,strcmp(uniqGenera,genera{i}))=blastOcc;
                        % firstOccurrenceMatrix(z1,z,j,k,i)=1;

                        subBigScores(z,j,k,i)=score;
                        subBigQCovS(z,j,k,i)=qCovS;
                        subBigPIdent(z,j,k,i)=pIdent;
                        subBigMatchNames{z,j,k,i}=matchName;
                        if subFirstOccurrenceMatrix(z,j,k,i) && blastOcc ~= subBlastOccurrenceMatrix(z,strcmp(uniqGenera,genera{i}))
                            disp(sprintf('WARNING DIFFERENT BLAST OCC: %d %d %d %d %d',z1,z,j,k,i));
                        end
                        subBlastOccurrenceMatrix(z,strcmp(uniqGenera,genera{i}))=blastOcc;
                        subFirstOccurrenceMatrix(z,j,k,i)=1;
                    end
                end
            end
        end
    end
    bigScores(z1,:,:,:,:) = subBigScores;
    bigQCovS(z1,:,:,:,:) = subBigQCovS;
    bigPIdent(z1,:,:,:,:) = subBigPIdent;
    blastOccurrenceMatrix(z1,:,:) = subBlastOccurrenceMatrix;
    firstOccurrenceMatrix(z1,:,:,:,:) = subFirstOccurrenceMatrix;
    bigMatchNames(z1,:,:,:,:) = subBigMatchNames;
end