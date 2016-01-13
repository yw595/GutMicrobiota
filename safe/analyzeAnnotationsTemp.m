configGutMicrobiota;
outputDir1 = [outputDir filesep 'analyzeAnnotations'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

avoidGenera={'Allobaculum','Clostridium','Butyricimonas', 'Oscillibacter','candidate_division_TM7','Rothia','Shigella'};
uniqGeneraTemp = uniqGenera(~ismember(uniqGenera,avoidGenera));

genesToLinesAnnot = containers.Map;
genesToLines = containers.Map;
for j=1:length(complexesAnnot)
    genes=complexesToGenesAnnot(complexesAnnot{j});
    for k=1:length(genes)
        FI=fopen([extraDir filesep complexesAnnot{j} '_' genes{k} '.txt']);
        lines = textscan(FI,'%s','Delimiter','\n');
        lines = lines{1};
        fclose(FI);
        genesToLinesAnnot([complexesAnnot{j} '_' genes{k}]) = lines;
        for l=1:length(lines)
            %lines{l}=regexptranslate('wildcard',regexptranslate('escape',lines{l}));
            lines{l} = strrep(lines{l},' ','_');
        end
        genesToLines([complexesAnnot{j} '_' genes{k}]) = lines;
    end
end

allGenes = keys(genesToLines);
occurrenceMatrix=[];%zeros(length(allGenera),length(allGenes));
for i=1:length(allGenera)
    disp(allGenera{i})
    %count=0;
    grepCmdPrefix = ['grep -P "\d\|('];
    for j=1:length(allGenes)
        lines = genesToLines(allGenes{j});
        grepCmd=grepCmdPrefix;
        for k=1:length(lines)-1
            grepCmd=[grepCmd lines{k} '|'];
            %occurrenceMatrix(i,j)=1;
        end
        grepCmd=[grepCmd lines{end} ')$" '];
        grepCmd=[grepCmd extraDir filesep allGenera{i} '.faa'];
        %disp(grepCmd);
        [status output]=system(grepCmd);
        %disp(output)
        if ~strcmp(output,'')
            disp(output)
            %disp(j)
            %disp(occurrenceMatrix)
            occurrenceMatrix(i,j)=1;
        end
    end
    %FI=fopen([extraDir filesep allGenera{i} '.faa']);
    %line = fgetl(FI);
    %while line~=-1
        % if regexp(line,'>')
        %     count=count+1;
        %     for j=1:length(allGenes)
        %         lines = genesToLines(allGenes{j});
        %         for k=1:length(lines)
        %             if ~isempty(regexp(line,lines{k}))
        %                 occurrenceMatrix(i,j)=1;
        %             end
        %         end
        %     end
        %     disp(count)
        % end
        %line=fgetl(FI);
    %end
    %fclose(FI);
end

%allGenes = keys(genesToLines);

% complexesToGenes=containers.Map({'NADH_dehydrogenase','cytochrome_bd_oxidase','cytochrome_bo_oxidase','cytochrome_c_oxidase','cytochrome_d_oxidase','cytochrome_o_oxidase','ATP_synthase'},{{'A','B','C','D','E','F','G','H','I','J','K','L','M','N'},{''},{''},{'I','II','III','IV'},{'I','II'},{'I','II','III','IV'},{'A','B','C','D','E'}});
% complexes = keys(complexesToGenes);

% extraDir = '/mnt/extra/blast/HMPRef';
% genesToLines = containers.Map;
% for j=1:length(complexes)
%     genes=complexesToGenes(complexes{j});
%     for k=1:length(genes)
%         if ~strcmp(genes{k},'')
%             FI=fopen([extraDir filesep complexes{j} '_' genes{k} '.txt']);
%         else
%             FI=fopen([extraDir filesep complexes{j} '.txt']);
%         end
%         lines = textscan(FI,'%s','Delimiter','\n');
%         lines = lines{1};
%         for l=1:length(lines)
%             lines{l}=regexptranslate('wildcard',regexptranslate('escape',lines{l}));
%         end
%         fclose(FI);
%         genesToLines([complexes{j} '_' genes{k}]) = lines;
%     end
% end
% allGenes = keys(genesToLines);

% occurrenceMatrix=[];%zeros(length(allGenera),length(allGenes));
% for i=1:length(allGenera)
%     disp(allGenera{i})
%     %count=0;
%     grepCmdPrefix = ['grep -P "\d ('];
%     for j=1:length(allGenes)
%         lines = genesToLines(allGenes{j});
%         grepCmd=grepCmdPrefix;
%         for k=1:length(lines)-1
%             grepCmd=[grepCmd lines{k} '|'];
%             %occurrenceMatrix(i,j)=1;
%         end
%         grepCmd=[grepCmd lines{end} ')$" '];
%         grepCmd=[grepCmd extraDir filesep allGenera{i} '.faa'];
%         %disp(grepCmd);
%         [status output]=system(grepCmd);
%         %disp(output)
%         if ~strcmp(output,'')
%             disp(output)
%             %disp(j)
%             %disp(occurrenceMatrix)
%             occurrenceMatrix(i,j)=1;
%         end
%     end
%     %FI=fopen([extraDir filesep allGenera{i} '.faa']);
%     %line = fgetl(FI);
%     %while line~=-1
%         % if regexp(line,'>')
%         %     count=count+1;
%         %     for j=1:length(allGenes)
%         %         lines = genesToLines(allGenes{j});
%         %         for k=1:length(lines)
%         %             if ~isempty(regexp(line,lines{k}))
%         %                 occurrenceMatrix(i,j)=1;
%         %             end
%         %         end
%         %     end
%         %     disp(count)
%         % end
%         %line=fgetl(FI);
%     %end
%     %fclose(FI);
% end

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