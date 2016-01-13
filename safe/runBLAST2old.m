configGutMicrobiota;
outputDir1 = [outputDir filesep 'runBLASTTemp'];

ZhangZhaoUpGenera = changeNames(ZhangZhaoUpGenera);
ZhangZhaoDownGenera = changeNames(ZhangZhaoDownGenera);
ForslundhildebrandUpGenera = changeNames(ForslundHildebrandUpGenera);
ForslundHildebrandDownGenera = changeNames(ForslundHildebrandDownGenera);

complexesToGenes = containers.Map({'ndh','Complex I','cytochrome bd-II oxidase','ATPase','B. sub cytochrome c oxidase','cytochrome bd-I oxidase','cytochrome bo oxidase','B. sub cytochrome b/c1 oxidase','Complex II'}, ...
    {{'ndh'},{'nuoA','nuoB','nuoC','nuoE','nuoF','nuoG','nuoH','nuoI','nuoJ','nuoK','nuoL','nuoM','nuoN'}, ...
    {'appB','appC'},{'atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH'},{'ctaB','ctaC','ctaD','ctaE','ctaF'}, ...
    {'cydA','cydB','cydX'}, {'cyoA','cyoB','cyoC','cyoD'}, {'qcrA','qcrB','qcrC'}, {'sdhA','sdhB','sdhC','sdhD'}});
complexes = keys(complexesToGenes);

ZhangZhaoComplexesToEVals1 = [];
ZhangZhaoComplexesToEVals2 = {};
ForslundHildebrandComplexesToEVals1 = [];
ForslundHildebrandComplexesToEVals2 = {};

for z1 = 1:length(complexes)
    geneNames = complexesToGenes(complexes{z1});
    for z=1:length(geneNames)
        outputDir4 = [outputDir1 filesep geneNames{z}];
        
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
                
                extraDir = '/mnt/extra/blast';
                namesDmp = [extraDir filesep 'names.dmp'];
                namesDmpFind = [outputDir3 filesep 'names.dmp.find'];
                
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
            if j==1
                %titleString = [titleString 'ZhangZhao'];
            else
                %titleString = [titleString 'ForslundHildebrand'];
            end
            
            if ~all(scores < 0) %&& strcmp(geneNames{z},'appB')
                
                upNum = scores(1:length(upGenera)); upNum = upNum(upNum ~= -2);
                downNum = scores(length(upGenera)+1:end); downNum = downNum(downNum ~= -2);
                upP = sum(upNum < 1)/length(upNum);
                downP = sum(downNum < 1)/length(downNum);
                overallP = (length(upNum)*upP+length(downNum)*downP)/(length(upNum)+length(downNum));
                zScore = (downP-upP)/sqrt(overallP*(1-overallP)*(1/length(upNum)+1/length(downNum)));
                pVal = 1-normcdf(zScore, 0, 1);
                pValHyge = 1-hygecdf(sum(downNum < 1), length(upNum)+length(downNum), sum(downNum < 1)+sum(upNum < 1), length(downNum));
                
                disp(titleString)
                disp(pVal)
                figure('Visible','off');
                dispYVals = yvals; dispYVals(dispYVals < 0) = 0;
                bar(xvals,dispYVals,'b','FaceColor','b');
                ylabel(ylabelString,'FontSize',20);
                ylim([-1*max(yvals(:)) 2.75*max(yvals(:))]);
                xlim([0 length(xvals)]);
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
                set(gca,'XTickLabel',{});
                clear line;
                line([0 length(xlabels)], [0 0]);
                line([length(upGenera)+.5 length(upGenera)+.5], [-1*max(yvals(:)) 1.3*max(yvals(:))], 'LineWidth', 2, 'Color', 'r');
                dim = [.4 .8 .07 .07];
                str = { ['Num BLAST Matches in Met/Ber-Up: ' num2str(sum(upNum < 1)) '/' num2str(length(upNum))], ...
                    ['Num BLAST Matches in Met/Ber-Down: ' num2str(sum(downNum < 1)) '/' num2str(length(downNum))], ...
                    ['Binomial Test P-Val for Different Means: ' num2str(pVal)], ...
                    ['Hypergeometric Test P-Val for Different Means: ' num2str(pValHyge)] };
                annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);
                if j==1
                    studyDir = 'ZhangZhao';
                else
                    studyDir = 'ForslundHildebrand';
                end
                if ~exist([outputDir1 filesep 'pngFiles' filesep studyDir],'dir')
                    mkdir([outputDir1 filesep 'pngFiles' filesep studyDir]);
                end
                if ~exist([outputDir1 filesep 'pngFiles' filesep studyDir filesep complexes{z1}],'dir')
                    mkdir([outputDir1 filesep 'pngFiles' filesep studyDir filesep complexes{z1}]);
                end
                saveas(gcf,[outputDir1 filesep 'pngFiles' filesep studyDir filesep complexes{z1} filesep titleString '.png']);
                close(gcf);
            end
        end
    end
end

figure('Visible','off');
boxplot(ZhangZhaoComplexesToEVals1, ZhangZhaoComplexesToEVals2);
prevLabels = unique(ZhangZhaoComplexesToEVals2);
set(gca,'XTickLabel',{' '});
ylim([-1.3*max(ZhangZhaoComplexesToEVals1) 1.3*max(ZhangZhaoComplexesToEVals1)]);
for i=1:length(prevLabels)
    text(i,-1.3*max(ZhangZhaoComplexesToEVals1),prevLabels{i},'Rotation',90,'FontSize',10);
    if mod(i,2)==0
        line([i+.5 i+.5],[-1.3*max(ZhangZhaoComplexesToEVals1) 1.3*max(ZhangZhaoComplexesToEVals1)]);
        [~, pValT] = ttest2(ZhangZhaoComplexesToEVals1(strcmp(ZhangZhaoComplexesToEVals2,prevLabels{i})), ...
            ZhangZhaoComplexesToEVals1(strcmp(ZhangZhaoComplexesToEVals2,prevLabels{i-1})));
        text(i-1.5,max(ZhangZhaoComplexesToEVals1)*1.2,'ttest p-val','FontSize',10);
        text(i-1.5,max(ZhangZhaoComplexesToEVals1)*1.1,'for diff. ','FontSize',10);
        text(i-1.5,max(ZhangZhaoComplexesToEVals1)*1.0,'means','FontSize',10);
        text(i-1.5,max(ZhangZhaoComplexesToEVals1)*0.9,num2str(pValT),'FontSize',10);
    end   
end
xlabel('Complex','FontSize',20)
ylabel('E-value','FontSize',20)
saveas(gcf,[outputDir1 filesep 'pngFiles' filesep 'ZhangZhaoEVals.png']);
close(gcf);

figure('Visible','off');
boxplot(ForslundHildebrandComplexesToEVals1, ForslundHildebrandComplexesToEVals2);
prevLabels = unique(ZhangZhaoComplexesToEVals2);
set(gca,'XTickLabel',{' '});
ylim([-1.3*max(ForslundHildebrandComplexesToEVals1) 1.3*max(ForslundHildebrandComplexesToEVals1)]);
for i=1:length(prevLabels)
    text(i,-1.3*max(ForslundHildebrandComplexesToEVals1(:)),prevLabels{i},'Rotation',90,'FontSize',10);
    if mod(i,2)==0
        line([i+.5 i+.5],[-1.3*max(ForslundHildebrandComplexesToEVals1(:)) 1.3*max(ForslundHildebrandComplexesToEVals1(:))]);
        [~, pValT] = ttest2(ForslundHildebrandComplexesToEVals1(strcmp(ForslundHildebrandComplexesToEVals2,prevLabels{i})), ...
            ForslundHildebrandComplexesToEVals1(strcmp(ForslundHildebrandComplexesToEVals2,prevLabels{i-1})));
        text(i-1.5,max(ForslundHildebrandComplexesToEVals1)*1.2,'ttest p-val','FontSize',10);
        text(i-1.5,max(ForslundHildebrandComplexesToEVals1)*1.1,'for diff. ','FontSize',10);
        text(i-1.5,max(ForslundHildebrandComplexesToEVals1)*1.0,'means','FontSize',10);
        text(i-1.5,max(ForslundHildebrandComplexesToEVals1)*0.9,num2str(pValT),'FontSize',10);
    end
end
xlabel('Complex','FontSize',20)
ylabel('E-value','FontSize',20)
saveas(gcf,[outputDir1 filesep 'pngFiles' filesep 'ForslundHildebrandEVals.png']);
close(gcf);