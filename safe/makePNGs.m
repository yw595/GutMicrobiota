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
elseif useHMPGenomes
    outputDir1 = [outputDir1 'HMPRef'];
end
outputDir1 = [outputDir1 filesep 'pngFiles'];
                  
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

ZhangZhaoUpGenera = changeNames(ZhangZhaoUpGenera);
ZhangZhaoDownGenera = changeNames(ZhangZhaoDownGenera);
ForslundhildebrandUpGenera = changeNames(ForslundHildebrandUpGenera);
ForslundHildebrandDownGenera = changeNames(ForslundHildebrandDownGenera);

if runCDD
    complexesToGenes = containers.Map;
    CDDDir = '/mnt/extra/blast/CDD';
    complexFiles = {'NADH_dehydrogenase_CDD.temp.txt', 'cytochrome_CDD.temp.txt','ATP_synthase_CDD.temp.txt'};
    for i=1:length(complexFiles)
        FI1 = fopen([CDDDir filesep complexFiles{i}]);
        complexName = complexFiles{i};
        complexName = complexName(1:regexp(complexName,'_CDD')-1);
        complexCdds = {};
        line = fgetl(FI1);
        while line~=-1
            words = strsplit(line,'\t');
            cddGi = words{1};
            cddID = words{2}
            cddName = words{4};
            cddName = strrep(cddName,' ','_');
            cddName = strrep(cddName,';','');
            cddName = strrep(cddName,'\.','');
            cddName = [cddID '_' cddName];
            if isempty(regexp(cddName,''''))
                complexCdds{end+1} = cddName;
            end
            line = fgetl(FI1);
        end
        complexesToGenes(complexName) = complexCdds;
        fclose(FI1);
    end
elseif runBLAST
    complexesToGenes = containers.Map({'mitochondrial_complex_I','ndh','Complex I','cytochrome bd-II oxidase','ATPase','B. sub cytochrome c oxidase','cytochrome bd-I oxidase','cytochrome bo oxidase','B. sub cytochrome b/c1 oxidase','Complex II'}, ...
                                  {{'NDUFS1','NDUFS2','NDUFS3','NDUFS7','NDUFS8', 'NDUFV1','NDUFV2','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'}, ...
                   {'ndh'},{'nuoA','nuoB','nuoC','nuoE','nuoF','nuoG','nuoH','nuoI','nuoJ','nuoK','nuoL','nuoM','nuoN'}, ...
    {'appB','appC'},{'atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH'},{'ctaB','ctaC','ctaD','ctaE','ctaF'}, ...
    {'cydA','cydB','cydX'}, {'cyoA','cyoB','cyoC','cyoD'}, {'qcrA','qcrB','qcrC'}, {'sdhA','sdhB','sdhC','sdhD'}});
end

complexes = keys(complexesToGenes);

ZhangZhaoComplexesToEVals1 = [];
ZhangZhaoComplexesToEVals2 = {};
ForslundHildebrandComplexesToEVals1 = [];
ForslundHildebrandComplexesToEVals2 = {};

for z1 = 1:length(complexes)
    geneNames = complexesToGenes(complexes{z1});
    for z=1:length(geneNames)
        if runCDD
            outputDir0 = [outputDir filesep 'runCDD'];
        elseif runBLAST
            outputDir0 = [outputDir filesep 'runBLAST'];
        end
        outputDir4 = [outputDir0 filesep geneNames{z}];
        
        for j=1:2
            if j==1
                upGenera = ZhangZhaoUpGenera;
                downGenera = ZhangZhaoDownGenera;
                outputDir2 = [outputDir4 filesep 'ZhangZhao'];
            else
                upGenera = ForslundHildebrandUpGenera;
                downGenera = ForslundHildebrandDownGenera;
                outputDir2 = [outputDir4 filesep 'ForslundHildebrand'];
            end
            
            scores = [];
            for k=1:2
                if k==1
                    generaToSearch = upGenera;
                    outputDir3 = [outputDir2 filesep 'up'];
                else
                    generaToSearch = downGenera;
                    outputDir3 = [outputDir2 filesep 'down'];
                end
                
                for i=1:length(generaToSearch)
                    FI = fopen([outputDir3 filesep generaToSearch{i} '.blast']);
                    if FI ~= -1
                        line = fgetl(FI);
                        if line ~= -1
                            words = strsplit(line,'\t');
                            scores(end+1) = str2num(words{1});
                            if j==1
                                if k==1
                                    ZhangZhaoComplexesToEVals1(end+1) = str2num(words{1});
                                    ZhangZhaoComplexesToEVals2{end+1} = [complexes{z1} ' UP'];
                                else
                                    ZhangZhaoComplexesToEVals1(end+1) = str2num(words{1});
                                    ZhangZhaoComplexesToEVals2{end+1} = [complexes{z1} ' DOWN'];
                                end
                            else
                                if k==1
                                    ForslundHildebrandComplexesToEVals1(end+1) = str2num(words{1});
                                    ForslundHildebrandComplexesToEVals2{end+1} = [complexes{z1} ' UP'];
                                else
                                    ForslundHildebrandComplexesToEVals1(end+1) = str2num(words{1});
                                    ForslundHildebrandComplexesToEVals2{end+1} = [complexes{z1} ' DOWN'];
                                end
                            end
                        else
                            scores(end+1) = -1;
                        end
                        fclose(FI);
                    else
                        scores(end+1) = -2;
                    end
                end
            end
            
            xvals = 1:length(scores); yvals = scores; titleString = geneNames{z}; xlabels = [upGenera; downGenera]; ylabelString = 'e-Value';
            
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
                %disp(upP)
                %disp(downP)
                %disp(overallP)
                %disp(upNum)
                %disp(downNum)
                %disp(scores)
                zScore = (downP-upP)/sqrt(overallP*(1-overallP)*(1/length(upNum)+1/length(downNum)));
                pVal = 1-normcdf(zScore, 0, 1);
                pValHyge = 1-hygecdf(sum(downNum < 1), length(upNum)+length(downNum), sum(downNum < 1)+sum(upNum < 1), length(downNum));
                
                disp(titleString)
                disp(pVal)
                figure('Visible','off');
                dispYVals = yvals; dispYVals(dispYVals < 0) = 0;
                bar(xvals,dispYVals,'b','FaceColor','b');
                ylabel(ylabelString,'FontSize',20);
                if max(yvals(:))==0 && min(yvals(:)==0)
                    ylim([-10 10]);
                else
                    ylim([-1*max(yvals(:)) 2.75*max(yvals(:))]);
                end
                xlim([0 length(xvals)]);
                set(gca,'XTickLabel',{});
                title(titleString,'FontSize',20);
                for k=1:length(xlabels)
                    text(k,-1*max(yvals(:)),xlabels{k},'Rotation',90,'FontSize',10);
                end
                negIdxs = find(yvals < 0);
                for k=1:length(negIdxs)
                    if yvals(negIdxs(k))==-2
                        text(negIdxs(k),0,'No reference sequences to BLAST against','Rotation',90,'FontSize',10);
                    elseif yvals(negIdxs(k))==-1
                        text(negIdxs(k),0,'No good BLAST matches','Rotation',90,'FontSize',10);
                    end
                end
                text(length(upGenera)/2,1.3*max(yvals(:)),'Met/Ber-Up','FontSize',10);
                text(length(upGenera)+length(downGenera)/2,1.3*max(yvals(:)),'Met/Ber-Down','FontSize',10);
                clear line;
                line([0 length(xlabels)], [0 0]);
                line([length(upGenera)+.5 length(upGenera)+.5], [-1*max(yvals(:)) 1.3*max(yvals(:))], 'LineWidth', 2, 'Color', 'r');
                dim = [.4 .8 .07 .07];
                str = { ['Num BLAST Matches in Met/Ber-Up: ' num2str(sum(upNum < 1)) '/' num2str(length(upNum))], ...
                    ['Num BLAST Matches in Met/Ber-Down: ' num2str(sum(downNum < 1)) '/' num2str(length(downNum))], ...
                    ['Binomial Test P-Val for Different Means: ' num2str(pVal)], ...
                    ['Hypergeometric Test P-Val for Different Means: ' num2str(pValHyge)] };
                annotation('textbox',dim,'String',str, ...
                           'FitBoxToText','on','FontSize',10);
                
                saveas(gcf,[outputDir3 filesep titleString '.png']);
                close(gcf);
            end
        end
    end
end

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