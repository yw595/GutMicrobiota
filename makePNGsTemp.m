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

complexes = keys(complexesToQueries);
complexesToComplexes = containers.Map({'Complex I','ATPase','B. sub cytochrome c oxidase','B. sub cytochrome b/c1 oxidase','Complex II'},{'NADH_dehydrogenase','ATP_synthase','cytochrome_c_oxidase','cytochrome_c_oxidase','succinate_dehydrogenase'});
matchingComplexes = keys(complexesToComplexes);

bigScores = -3*ones(length(complexes),50,2,4,20);
bigMatchNames = {};
eValThresh = 0.00001;
for z1 = 1:length(complexes)
    geneNames = complexesToQueries(complexes{z1});
    for z=1:length(geneNames)
        outputDir4 = [outputDir0 filesep geneNames{z}];
        
        for j=1:2
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
            
            scores = [];
            for k=1:4
                if k==1
                    generaToSearch = upGenera;
                    outputDir3 = [outputDir2 filesep 'up'];
                elseif k==2
                    generaToSearch = downGenera;
                    outputDir3 = [outputDir2 filesep 'down'];
                elseif k==3
                    generaToSearch = upStrictGenera;
                    outputDir3 = [outputDir2 filesep 'upStrict'];
                else
                    generaToSearch = downStrictGenera;
                    outputDir3 = [outputDir2 filesep 'downStrict'];
                end
                
                for i=1:length(generaToSearch)
                    FI = fopen([outputDir3 filesep generaToSearch{i} '.blast']);
                    if FI ~= -1
                        line = fgetl(FI);
                        if line ~= -1
                            words = strsplit(line,'\t');
                            scores(end+1) = str2num(words{1});
                            bigScores(z1,z,j,k,i) = str2num(words{1});
                            bigMatchNames{z1,z,j,k,i} = words{9};
                        else
                            scores(end+1) = -1;
                            bigScores(z1,z,j,k,i) = -1;
                            bigMatchNames{z1,z,j,k,i} = 'NA';
                        end
                        fclose(FI);
                    else
                        scores(end+1) = -2;
                        bigScores(z1,z,j,k,i) = -2;
                        bigMatchNames{z1,z,j,k,i} = 'NA';
                    end
                end
            end
            if 0
            xvals = 1:length(scores); yvals = scores; titleString = geneNames{z}; xlabels = [upGenera downGenera]; ylabelString = 'e-Value'; 
            labels = {};
            for tempIter=1:length(upGenera)
                labels{tempIter}='Met/Ber Up';
            end
            for tempIter=1:length(downGenera)
                labels{end+1}='Met/Ber Down';
            end
            
            if ~all(scores < 0) %&& strcmp(geneNames{z},'NDUFS8')
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
                
                upNum = scores(1:length(upGenera)); upNum = upNum(upNum ~= -2);
                downNum = scores(length(upGenera)+1:end); downNum = downNum(downNum ~= -2);
                upP = sum(upNum < 1)/length(upNum);
                downP = sum(downNum < 1)/length(downNum);
                overallP = (length(upNum)*upP+length(downNum)* downP)/(length(upNum)+length(downNum));
                zScore = (downP-upP)/sqrt(overallP*(1-overallP)*(1/length(upNum)+1/length(downNum)));
                pVal = 1-normcdf(zScore, 0, 1);
                pValHyge = 1-hygecdf(sum(downNum < 1), length(upNum)+length(downNum), sum(downNum < 1)+sum(upNum < 1), length(downNum));
                
                disp(titleString)
                disp(pVal)

                specialVals = containers.Map([-2 -1],{'No reference sequences to BLAST against','No good BLAST matches'});
                
                % dim = [.4 .8 .07 .07];
                % str = { ['Num BLAST Matches in Met/Ber-Up: ' num2str(sum(upNum < 1)) '/' num2str(length(upNum))], ...
                %     ['Num BLAST Matches in Met/Ber-Down: ' num2str(sum(downNum < 1)) '/' num2str(length(downNum))], ...
                %     ['Binomial Test P-Val for Different Means: ' num2str(pVal)], ...
                %     ['Hypergeometric Test P-Val for Different Means: ' num2str(pValHyge)] };
                % annotation('textbox',dim,'String',str, ...
                %            'FitBoxToText','on','FontSize',10);
                makeBarGut(xvals,yvals,titleString,xlabels,ylabelString,outputDir3,labels,specialVals);
            end
            end
        end
    end
end

if 1
generaMatrix = {'ZhangZhaoUpGenera' 'ZhangZhaoDownGenera' 'ZhangZhaoUpStrictGenera' 'ZhangZhaoDownStrictGenera'; 'ForslundHildebrandUpGenera' 'ForslundHildebrandDownGenera' 'ForslundHildebrandUpStrictGenera' 'ForslundHildebrandDownStrictGenera'};
for studyIdx = 1:2
    for upDownIdx = 1:4
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
            %scores = [];
            xvals = [];
            yvals = [];
            xlabels = {};
            labels = {};
            blastPers = zeros(length(matchingComplexes),1);
            annotPers = blastPers;
            complexThere = zeros(length(complexes),1);
            for l=1:length(complexes)
                geneNames = complexesToQueries(complexes{l});
                genesThere = zeros(length(geneNames),1);
                for m=1:length(geneNames)
                    %scores(l,m) =
                    %bigScores(l,m,studyIdx,upDownIdx, k);
                    xvals(end+1) = length(xvals)+1;
                    yvals(end+1) = bigScores(l,m,studyIdx,upDownIdx, k);
                    xlabels{end+1} = sprintf([geneNames{m} '\t' bigMatchNames{l,m,studyIdx,upDownIdx,k}]);
                    labels{end+1} = complexes{l};
                    genesThere(m) = bigScores(l,m,studyIdx,upDownIdx, k)<=eValThresh;
                end
                complexThere(k) = genesThere
                
                if any(strcmp(complexes{l},matchingComplexes)) && any(strcmp(allGenera,genera{k}))
                    shortGeneNames = complexesToGenesAnnot(complexesToComplexes(complexes{l}));
                    longGeneNames = {};
                    for tempIdx=1:length(shortGeneNames)
                        longGeneNames{tempIdx} = [complexesToComplexes(complexes{l}) '_' shortGeneNames{tempIdx}];
                    end
                    blastPers(strcmp(complexes{l},matchingComplexes)) = sum(yvals<1 & yvals>=0)/length(yvals);
                    occurrences = occurrenceMatrix(strcmp(allGenera,genera{k}),ismember(allGenes,longGeneNames));
                    annotPers(strcmp(complexes{l},matchingComplexes)) = sum(occurrences)/length(occurrences);
                end
            end
            
            titleString = genera{k}; ylabelString = 'e-Value';            
            specialVals = containers.Map([-2 -1],{'No reference sequences to BLAST against','No good BLAST matches'});            
            makeBarGut(xvals,yvals,titleString,xlabels,ylabelString,outputDirTemp,labels,specialVals);

            makeBarGut(1:length(annotPers),[blastPers annotPers], [titleString 'Comp'],matchingComplexes,'Percentage',outputDirTemp,{},specialVals);
        end
    end
end
end

if 0
for j=1:2
    if j==1
        complexesToEVals1 = ZhangZhaoComplexesToEVals1;
        complexesToEVals2 = ZhangZhaoComplexesToEVals2;
    else
        complexesToEVals1 = ForslundHildebrandComplexesToEVals1;
        complexesToEVals2 = ForslundHildebrandComplexesToEVals2;
    end
    figure('Visible','off');
    boxplot(complexesToEVals1, complexesToEVals2);
    prevLabels = unique(complexesToEVals2);
    set(gca,'XTickLabel',{' '});
    ylim([-1.3*max(complexesToEVals1) 1.3*max(complexesToEVals1)]);
    for i=1:length(prevLabels)
        text(i,-1.3*max(complexesToEVals1),prevLabels{i},'Rotation',90,'FontSize',10);
        if mod(i,2)==0
            clear line
            line([i+.5 i+.5],[-1.3*max(complexesToEVals1) 1.3*max(complexesToEVals1)]);
            [~, pValT] = ttest2(complexesToEVals1(strcmp(complexesToEVals2,prevLabels{i})), ...
                complexesToEVals1(strcmp(complexesToEVals2,prevLabels{i-1})));
            text(i-1.5,max(complexesToEVals1)*1.2,'ttest p-val','FontSize',10);
            text(i-1.5,max(complexesToEVals1)*1.1,'for diff. ','FontSize',10);
            text(i-1.5,max(complexesToEVals1)*1.0,'means','FontSize',10);
            text(i-1.5,max(complexesToEVals1)*0.9,num2str(pValT),'FontSize',10);
        end   
    end
    xlabel('Complex','FontSize',20)
    ylabel('E-value','FontSize',20)
    if j==1
        saveas(gcf,[outputDir1 filesep 'ZhangZhaoEVals.png']);
    else
        saveas(gcf,[outputDir1 filesep 'ForslundHildebrandEVals.png']);
    end
    close(gcf);
end
end