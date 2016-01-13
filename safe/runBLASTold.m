%speciesToSearch = {'E. coli O157:H7 Bacteriophage HY01','Parasutterella','Bilophila','Bacteroides', ...
%    'Alistipes','Blautia','Lachnospiraceae','Clostridium','Allobaculum','Holdemania','Butyricicoccus','Phascolarctobacterium'};%'Butyricimonas'};
actualRun = 1;
outputDir = '/home/ubuntu/MATLAB/GutMicrobiota/output';
outputDir1 = [outputDir filesep 'runBLAST'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

geneNames = {'nuoA','nuoB','nuoC','nuoE','nuoF','nuoG','nuoH','nuoI','nuoJ','nuoK','nuoL','nuoM','nuoN', ...
    'appB','appC','atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH','ctaB','ctaC','ctaD','ctaE','ctaF', ...
    'cydA','cydB','cydX','cyoA','cyoB','cyoC','cyoD','qcrA','qcrB','qcrC','sdhA','sdhB','sdhC','sdhD'};

for z=1:length(geneNames)
    outputDir4 = [outputDir1 filesep geneNames{z}];
    if ~exist(outputDir4,'dir')
        mkdir(outputDir4)
    end
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
        
        if ~exist(outputDir2,'dir')
            mkdir(outputDir2)
        end
        
        for l=1:2
            if l==1
                generaTemp = upGenera;
            else
                generaTemp = downGenera;
            end
            for k=1:length(generaTemp)
                if ~isempty(regexp(generaTemp{k},'Clostridium'))
                    generaTemp{k} = 'Clostridium';
                end
                if ~isempty(regexp(generaTemp{k},'Escherichia/Shigella'))
                    generaTemp{k} = 'Escherichia';
                end
                if ~isempty(regexp(generaTemp{k},'TM7'))
                    generaTemp{k} = 'candidate division TM7';
                end
                if ~isempty(regexp(generaTemp{k},'Lachnospiracea'))
                    generaTemp{k} = 'Lachnospiraceae';
                end
            end
            if l==1
                upGenera = generaTemp;
            else
                downGenera = generaTemp;
            end
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
            
            if ~exist(outputDir3,'dir')
                mkdir(outputDir3)
            end
            
            %outputDir = '/mnt/extra/blast';
            extraDir = '/mnt/extra/blast';
            namesDmp = [extraDir filesep 'names.dmp'];
            namesDmpFind = [outputDir3 filesep 'names.dmp.find'];
            command = 'grep -P "\t(';
            for i=1:length(generaToSearch)-1
                command = [command generaToSearch{i} '|'];
            end
            command = [command generaToSearch{end} ')\t" ' namesDmp ' | cut -f3,1 > ' namesDmpFind];
            disp(command)
            % do not check actualRun here, because this is needed for reading
            % namesDmpFind below
            status=system(command);
            
            FI = fopen(namesDmpFind);
            line = fgetl(FI);
            namesToTaxids = containers.Map;
            while line ~= -1
                words = strsplit(line,'\t');
                namesToTaxids(words{2}) = words{1};
                line = fgetl(FI);
            end
            fclose(FI);
            
            command = sprintf(['gawk ''{if($2==%s) {print $1 > "' outputDir3 filesep '%s.txt"}'],namesToTaxids(generaToSearch{1}),generaToSearch{1});
            for i=2:length(generaToSearch)
                command = [command sprintf([' else if($2==%s) {print $1 > "' outputDir3 filesep '%s.txt"}'],namesToTaxids(generaToSearch{i}),generaToSearch{i})];
            end
            command = [command '}'' ' extraDir filesep 'gi_taxid_prot.dmp'];
            disp(command)
            command = ['gawk ''BEGIN {FS=OFS="\t"}; FNR==NR{array[$1]=$2; next}{if($2 in array){print $1; print "' outputDir3 filesep '" array[$2] ".txt"; print $1 > "' outputDir3 filesep '" array[$2] ".txt"}}'' ' outputDir3 filesep 'names.dmp.find ' extraDir filesep 'gi_taxid_prot.dmp'];
            disp(command);
            if actualRun
                status=system(command);
            end
            
            commandPrefix = ['blastp -db nr -query /home/ubuntu/MATLAB/GutMicrobiota/input' filesep geneNames{z} '.fasta -gilist "' outputDir3 filesep];
            for i=1:length(generaToSearch)
                command = [commandPrefix generaToSearch{i} '.txt" -outfmt "6 evalue qacc sacc qstart qend sstart send" -out "' outputDir3 filesep generaToSearch{i} '.blast"'];
                disp(command);
                if actualRun
                    status=system(command);
                end
            end
            
            for i=1:length(generaToSearch)
                FI = fopen([outputDir3 filesep generaToSearch{i} '.blast']);
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
        
        %makeBar(1:length(scores), scores, 'blastVals', generaToSearch, 'e-Value', outputDir);
        
        xvals = 1:length(scores); yvals = scores; titleString = 'blastVals'; xlabels = [upGenera; downGenera]; ylabelString = 'e-Value';
        if j==1
            titleString = [titleString 'ZhangZhao'];
        else
            titleString = [titleString 'ForslundHildebrand'];
        end
        
        figure('Visible','off');
        bar(xvals,yvals,'b','FaceColor','b');
        ylabel(ylabelString,'FontSize',20);
        ylim([-.1*max(yvals(:)) 1.3*max(yvals(:))]);
        title(titleString,'FontSize',20);
        for k=1:length(xlabels)
            text(k,-.1*max(yvals(:)),xlabels{k},'Rotation',90,'FontSize',20);
        end
        set(gca,'XTickLabel',{});
        clear line;
        line([0 length(xlabels)], [0 0]);
        line([length(upGenera) length(upGenera)], [-.1*max(yvals(:)) 1.3*max(yvals(:))]);
        saveas(gcf,[outputDir2 filesep titleString '.png']);
        close(gcf);
    end
end