configGutMicrobiota;
outputDir1 = [outputDir filesep 'analyzeAnnotations'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

genesToLinesAnnot = containers.Map;
genesToLines = containers.Map;
for j=1:length(complexesAnnot)
    genes=complexesToGenesAnnot(complexesAnnot{j});
    for k=1:length(genes)
        FI=fopen([HMPDir filesep complexesAnnot{j} '_' genes{k} '.txt']);
        lines = textscan(FI,'%s','Delimiter','\n');
        lines = lines{1};
        fclose(FI);
        for l=1:length(lines)
            lines{l} = strrep(lines{l},' ','_');
        end
        disp([complexesAnnot{j} '_' genes{k}])
        genesToLines(genesToGenes([complexesAnnot{j} '_' genes{k}])) = lines;
        for l=1:length(lines)
            lines{l}=regexptranslate('wildcard',regexptranslate('escape',lines{l}));
        end
        genesToLinesAnnot(genesToGenes([complexesAnnot{j} '_' genes{k}])) = lines;
    end
end

occurrenceMatrix=[];
for i=1:length(uniqGeneraAnnot)
    disp(uniqGeneraAnnot{i})
    grepCmdPrefix = ['grep -P "\d\|('];
    for j=1:length(matchingGenesArray)
        lines = genesToLinesAnnot(matchingGenesArray{j});
        grepCmd = grepCmdPrefix;
        for k=1:length(lines)-1
            grepCmd = [grepCmd lines{k} '|'];
        end
        grepCmd = [grepCmd lines{end} ')$" '];
        grepCmd = [grepCmd HMPDir filesep uniqGeneraAnnot{i} '.faa'];
        [status output]=system(grepCmd);
        if ~strcmp(output,'')
            disp(output)
            occurrenceMatrix(i,j)=1;
        end
    end
end

% outputDir1 = [outputDir1 filesep 'pngFilesNames'];
% if ~exist(outputDir1,'dir')
%     mkdir(outputDir1)
% end
% for i=1:length(allGenera)
%     figure('Visible','off');
%     bar(1:length(allGenes),occurrenceMatrix(i,:),'b','FaceColor', 'b');
%     ylim([-10,2])
%     for j=1:length(allGenes)
%         text(j,-10,allGenes{j},'FontSize',10,'Rotation',90);
%     end
%     title(allGenera{i});
%     saveas(gcf,[outputDir1 filesep allGenera{i} '.png']);
%     close(gcf);
% end