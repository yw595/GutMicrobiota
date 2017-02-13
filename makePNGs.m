configGutMicrobiota;
runCDD = 0;
runBLAST = 1;
if runCDD
    outputDir1 = [outputDir filesep 'runCDD'];
elseif runBLAST
    outputDir1 = [outputDir filesep 'runBLAST'];
end
useTempGenomes = 0;
useHMPGenomes = 1;
if useTempGenomes
    outputDir1 = [outputDir1 'Temp'];
elseif useHMPGenomes
    outputDir1 = [outputDir1 'HMP'];
end
outputDir0 = outputDir1;
outputDir1 = [outputDir1 filesep 'pngFiles'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

if runCDD
    complexesToQueries = complexesToCdds;
elseif runBLAST
    complexesToQueries = complexesToGenes;
end

%complexes = keys(complexesToQueries);
%matchingComplexes = keys(complexesToComplexes);

%otherGenera = keys(generaToFileNames);
%otherGenera = otherGenera(~ismember(otherGenera,uniqGenera));

eValThresh = 0.00001;
bigScoresOther = -3*ones(length(complexes),50,20);
bigQCovSOther = -3*ones(length(complexes),50,20);
bigPIdentOther = -3*ones(length(complexes),50,20);
for z1 = 1:length(complexes)
    geneNames = complexesToQueries(complexes{z1});
    for z=1:length(geneNames)
        outputDir4 = [outputDir0 filesep geneNames{z}];
        outputDir2 = [outputDir4 filesep 'other'];

        generaToSearch = otherGenera;
        %outputDir3 = [outputDir2 filesep 'other'];

        for i=1:length(generaToSearch)
            BLASTFile=[outputDir2 filesep generaToSearch{i} '.blast'];
            [score matchName qCovS pIdent blastOcc] = readBLASTFile(BLASTFile,0);
            bigScoresOther(z1,z,i)=score;
            bigQCovSOther(z1,z,i)=qCovS;
            bigPIdentOther(z1,z,i)=pIdent;
            %blastOccurrenceMatrixOther
            % if FI ~= -1
            %     line = fgetl(FI);
            %     if line ~= -1
            %         words = strsplitYiping(line,'\t');
            %         bigScoresOther(z1,z,i) = str2num(words{1});
            %     else
            %         bigScoresOther(z1,z,i) = -1;
            %     end
            %     fclose(FI);
            % else
            %     bigScoresOther(z1,z,i) = -2;
            % end
        end
    end
end
            
for z1 = 13:13%length(complexes)
    geneNames = complexesToQueries(complexes{z1});
    for z=1:1%length(geneNames)
        outputDir4 = [outputDir0 filesep geneNames{z}];
        
        for j=1:1%2
            if j==1
                upGenera = ZhangZhaoUpGenera;
                downGenera = ZhangZhaoDownGenera;
                upStrictGenera = ZhangZhaoUpStrictGenera;
                downStrictGenera = ZhangZhaoDownStrictGenera;
                outputDir2 = [outputDir4 filesep 'ZhangZhao'];
            else
                upGenera = ForslundHildebrandUpGenera;
                downGenera = ForslundHildebrandDownGenera;
                upStrictGenera = ForslundHildebrandUpStrictGenera;
                downStrictGenera = ForslundHildebrandDownStrictGenera;
                outputDir2 = [outputDir4 filesep 'ForslundHildebrand'];
            end
            if 1
            for scoresIdx = 1:4
                titleString = geneNames{z}    
                if scoresIdx==1
                    scoresMatrix=bigScores;specificThresh=0.00005;ylabelString = 'e-Value';titleString=[titleString ' E-Value'];scoresMatrixOther=bigScoresOther;
                elseif scoresIdx==2
                    scoresMatrix=bigQCovS;specificThresh=.3;ylabelString = 'Query-Coverage';titleString=[titleString ' Query-Coverage'];scoresMatrixOther=bigQCovSOther;
                elseif scoresIdx==3
                    scoresMatrix=bigPIdent;specificThresh=.3;ylabelString = 'Percent-Identity';titleString=[titleString ' Percent-Identity'];scoresMatrixOther=bigPIdentOther;
                else
                    scoresMatrix=blastOccurrenceMatrix;specificThresh=1;ylabelString = 'BLAST-Hit';titleString=[titleString ' BLAST-Hits'];scoresMatrixOther=bigScoresOther;
                end
                if scoresIdx~=4
                    scores = [squeeze(scoresMatrix(z1,z,j,3,1:length(upStrictGenera)))' squeeze(scoresMatrix(z1,z,j,4,1:length(downStrictGenera)))'];
                else
                    scores = [squeeze(scoresMatrix(z1,z,ismember(uniqGenera,upStrictGenera)))' squeeze(scoresMatrix(z1,z,ismember(uniqGenera,downStrictGenera)))'];
                end

                disp(scores)
                
                xvals = 1:length(scores); yvals = scores; xlabels = [upStrictGenera downStrictGenera]; 

                if ~all(scores < 0)
                    if j==1
                        outputDir2 = [outputDir1 filesep 'ZhangZhao'];
                    else
                        outputDir2 = [outputDir1 filesep 'ForslundHildebrand'];
                    end
                    if ~exist(outputDir2,'dir')
                        mkdir(outputDir2);
                    end
                    outputDir3 = [outputDir2 filesep complexes{z1}];
                    if ~exist(outputDir3,'dir')
                        mkdir(outputDir3);
                    end

                    upNum = scores(1:length(upStrictGenera)); upNum = upNum(upNum ~= -2);
                    downNum = scores(length(upStrictGenera)+1:end); downNum = downNum(downNum ~= -2);
                    if scoresIdx==2 || scoresIdx==3
                        upValid = sum(upNum > specificThresh & upNum >= 0);
                        downValid = sum(downNum > specificThresh & downNum >= 0);
                    else
                        upValid = sum(upNum < specificThresh & upNum >= 0);
                        downValid = sum(downNum < specificThresh & downNum >= 0);
                    end
                    upP = upValid/length(upNum);
                    downP = downValid/length(downNum);
                    overallP = (length(upNum)*upP+length(downNum)* downP)/(length(upNum)+length(downNum));
                    zScore = (downP-upP)/sqrt(overallP*(1-overallP)*(1/length(upNum)+1/length(downNum)));
                    pVal = 1-normcdf(zScore, 0, 1);
                    pValHyge = 1-hygecdf(downValid, length(upNum)+length(downNum), downValid+upValid, length(downNum));

                    disp(titleString)
                    disp(pVal)

                    specialVals = containers.Map([-2 -1],{'No reference sequences to BLAST against','No good BLAST matches'});

                    str = [sprintf('Num BLAST Matches in Met/Ber-Up: %d/%d ',sum(upNum < 1),length(upNum)), sprintf('Num BLAST Matches in Met/Ber-Down: %d\%d ',sum(downNum < 1),length(downNum)), sprintf('Binomial Test P-Val for Different Means: %f ',pVal), sprintf('Hypergeometric Test P-Val for Different Means: %f',pValHyge)];
                    labels = {};
                    for tempIter=1:length(upStrictGenera)
                        labels{tempIter}=['Met/Ber Up ' str];
                    end
                    for tempIter=1:length(downStrictGenera)
                        labels{end+1}=['Met/Ber Down ' str];
                    end
                    makeBar(xvals,yvals,titleString,outputDir3,'ylabelString',ylabelString,'xlabels',xlabels,'labels',labels,'specialVals',specialVals);

                    xvals = 1:100; yvals1 = yvals; yvals = squeeze(scoresMatrixOther(z1,z,:)); yvals = yvals(yvals>=0);
                    bins = linspace(min(yvals),max(yvals),100);
                    yvalsHist = histc(yvals,bins)';  
                    yvals1Hist = histc(yvals1,bins);
                    titleString = [titleString 'Hist']; labels = '';
                    if ~isempty(regexp(titleString,'Percent'))
                        disp('scoresMatrixOther')
                        disp(scoresMatrixOther(z1,z,:));
                        disp('bins')
                        disp(bins);
                        disp('yvals1')
                        disp(yvals1)
                        disp('yvals1Hist')
                        disp(yvals1Hist);
                        disp('xvals')
                        disp(xvals)
                        disp('yvalsHist')
                        disp(yvalsHist)
                        %absurd = absurd+1;
                    end
                    %yvalsHist(1)=0;
                    %yvalsHist(2)=yvalsHist(1);
                    %xvals = xvals(2:end); yvalsHist = yvalsHist(2:end);
                    makeBar(xvals,yvalsHist,titleString,outputDir3,'ylabelString','Frequency','indLines',yvals1Hist);
                end
            end
            end
        end
    end
end

complexCorrMatrix = zeros(length(complexes),length(complexes));
for z1 = 1:length(complexes)
    geneNames = complexesToQueries(complexes{z1});
    for i=1:length(uniqGenera)
        pcts(z1,i) = sum(blastOccurrenceMatrix(z1,1:length(geneNames),i))/length(blastOccurrenceMatrix(z1,1:length(geneNames),i));    
    end
end
for z1 = 1:length(complexes)
    for z2 = 1:length(complexes)
        [rho pval] = corr(pcts(z1,:)',pcts(z2,:)','type','Spearman');
        complexCorrMatrixRho(z1,z2) = rho;
        complexCorrMatrixPVal(z1,z2) = pval;
    end
end
for z1 = 1:length(complexes)
    for z2=z1:length(complexes)
        complexCorrMatrixRho(z1,z1)=0;
        complexCorrMatrixPVal(z1,z2)=Inf;
    end
end
[~, sortIdxs] = sort(complexCorrMatrixPVal(:),'descend');
for i=0:9
    [complexIdx1 complexIdx2] = ind2sub(size(complexCorrMatrixPVal), sortIdxs(end-i));
    disp([complexes{complexIdx1} ' ' complexes{complexIdx2}])
    xvals=pcts(complexIdx1,:);yvals=pcts(complexIdx2,:);
    titleString=sprintf('complexCorr %s vs %s R^2=%f p-val=%f',complexes{complexIdx1},complexes{complexIdx2},complexCorrMatrixRho(complexIdx1,complexIdx2),complexCorrMatrixPVal(complexIdx1,complexIdx2));
    ylabelString=[complexes{complexIdx1} 'Percent Coverage'];
    xlabelString=[complexes{complexIdx2} 'Percent Coverage'];
    if i==0
        %complexIdx1
        %complexIdx2
        %xvals
        %yvals
    end
    makeBar(xvals,yvals,titleString,outputDir1,'ylabelString',ylabelString,'isScatter',1,'xlabelString',xlabelString);
end

upStrictGeneraAll = unique([ZhangZhaoUpStrictGenera ForslundHildebrandUpStrictGenera]);
downStrictGeneraAll = unique([ZhangZhaoDownStrictGenera ForslundHildebrandDownStrictGenera]);
titleString = ['Complex_Diff_Coverage'];
xvals = []; yvals = [];
for i=1:length(complexes)
    pctsUp = pcts(i,ismember(uniqGenera,upStrictGeneraAll));
    pctsDown = pcts(i,ismember(uniqGenera,downStrictGeneraAll));
    xvals(i) = i; yvals(i,1) = mean(pctsUp); yvals(i,2) = mean(pctsDown);
    xlabels{i} = complexes{i}; [h p]=ttest2(pctsUp,pctsDown);
    labels{i} = ['ttest p-val: ' num2str(p)];
end
makeBar(xvals,yvals,titleString,outputDir1,'ylabelString','Average Coverage','xlabels',xlabels,'labels',labels,'xlabelsFontSize',25);

if 1
generaMatrix = {'ZhangZhaoUpGenera' 'ZhangZhaoDownGenera' 'ZhangZhaoUpStrictGenera' 'ZhangZhaoDownStrictGenera'; 'ForslundHildebrandUpGenera' 'ForslundHildebrandDownGenera' 'ForslundHildebrandUpStrictGenera' 'ForslundHildebrandDownStrictGenera'};
for studyIdx = 1:1
    for upDownIdx = 3:3
        if studyIdx==1
            outputDirTemp = [outputDir1 filesep 'ZhangZhao'];
        else
            outputDirTemp = [outputDir1 filesep 'ForslundHildebrand'];
        end
        if upDownIdx==1
            outputDirTemp = [outputDirTemp filesep 'up'];
        elseif upDownIdx==2
            outputDirTemp = [outputDirTemp filesep 'down'];
        elseif upDownIdx==3
            outputDirTemp = [outputDirTemp filesep 'upStrict'];    
        else
            outputDirTemp = [outputDirTemp filesep 'downStrict'];
        end
        if ~exist(outputDirTemp,'dir')
            mkdir(outputDirTemp);
        end
        
        genera = eval(generaMatrix{studyIdx,upDownIdx});
        for k=1:length(genera)
            xvals = [];
            yvals = [];
            xlabels = {};
            labels = {};
            blastPers = zeros(length(matchingComplexes),1);
            annotPers = blastPers;
            for l=1:length(complexes)
                geneNames = complexesToQueries(complexes{l});
                matchvals = [];
                for m=1:length(geneNames)
                    xvals(end+1) = length(xvals)+1;
                    yvals(end+1) = bigScores(l,m,studyIdx,upDownIdx, k);
                    matchvals(m) = blastOccurrenceMatrix(l,m,strcmp(uniqGenera,genera{k}));
                    %xlabels{end+1} = sprintf([geneNames{m} '\t' bigMatchNames{l,m,studyIdx,upDownIdx,k}]);
                    xlabels{end+1} = geneNames{m};
                    labels{end+1} = complexes{l};    
                end

                if any(strcmp(complexes{l},matchingComplexes)) && any(strcmp(uniqGeneraAnnot,genera{k}))
                    blastPers(strcmp(complexes{l},matchingComplexes)) = sum(matchvals)/length(matchvals);
                    occurrences = occurrenceMatrix(strcmp(uniqGeneraAnnot,genera{k}),ismember(matchingGenesArray,geneNames));
                    annotPers(strcmp(complexes{l},matchingComplexes)) = sum(occurrences)/length(occurrences);
                end
            end

            if strcmp(genera{k},'Escherichia')
                %disp(studyIdx)
                %disp(upDownIdx)
                %disp(annotPers)
                %disp(blastPers)
            end
            
            titleString = genera{k}; ylabelString = 'e-Value';            
            specialVals = containers.Map([-2 -1],{'No reference sequences to BLAST against','No good BLAST matches'});            
            makeBar(xvals,yvals,titleString,outputDirTemp,'ylabelString',ylabelString,'xlabels',xlabels,'labels',labels,'specialVals',specialVals);
            midpoint = floor(length(xvals)/2);
            xvals1 = xvals(1:midpoint); yvals1 = yvals(1:midpoint);
            titleString1 = [titleString '_1'];
            xlabels1 = xlabels(1:midpoint); labels1 = labels(1:midpoint);
            makeBar(xvals1,yvals1,titleString1,outputDirTemp,'ylabelString',ylabelString,'xlabels',xlabels1,'labels',labels1,'specialVals',specialVals);
            xvals2 = xvals(midpoint+1:end); yvals2 = yvals(midpoint+1:end);
            titleString2 = [titleString '_2'];
            xlabels2 = xlabels(midpoint+1:end); labels2 = labels(midpoint+1:end);
            makeBar(xvals2,yvals2,titleString2,outputDirTemp,'ylabelString',ylabelString,'xlabels',xlabels2,'labels',labels2,'specialVals',specialVals);
            
            makeBar(1:length(annotPers),[blastPers annotPers],[titleString 'Comp'],outputDirTemp,'ylabelString','Percentage','xlabels',matchingComplexes,'specialVals',specialVals);
        end
    end
end
end

if 0
for j=1:2
    if j==1
        titleString = 'ZhangZhaoEVals';
    else
        titleString = 'ForslundHildebrandEVals';
    end
    xvals = 1:length(complexes)*2;
    yvals = []; xlabels = {}; labels = {};
    for z1=1:length(complexes)
        geneNames = complexesToGenes(complexes{z1});
        for k=3:4
            genera = eval(generaMatrix{j,k});
            yvals = [yvals reshape(bigScores(z1,1:length(geneNames),j,k,1:length(genera)),1,length(geneNames)*length(genera))];
            if k==1
                prevyvals=yvals;
            end
            if k==1
                xlabels{end+1} = [complexes{z1} ' UP STRICT'];
            else
                xlabels{end+1} = [complexes{z1} ' DOWN STRICT'];
            end
            if k==2
                [~, pValT] = ttest2(prevyvals,yvals);
                newLabel = sprintf(['ttest p-val for diff. means: %f'],pValT);
                labels{end+1}=newLabel; labels{end+1}=newLabel;
            end
        end
    end
    ylabelString = 'e-Value';
    %makeBar(xvals,yvals,titleString,outputDir1,'xlabels',xlabels,'ylabelString',ylabelString,'labels',labels);
end
end