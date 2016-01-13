configGutMicrobiota;
outputDir1 = [outputDir filesep 'runBLAST'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

ZhangZhaoUpGenera = changeNames(ZhangZhaoUpGenera);
ZhangZhaoDownGenera = changeNames(ZhangZhaoDownGenera);
ForslundHildebrandUpGenera = changeNames(ForslundHildebrandUpGenera);
ForslundHildebrandDownGenera = changeNames(ForslundHildebrandDownGenera);

%geneNames = {'NDUFS1','NDUFS2','NDUFS3','NDUFS7','NDUFS8','NDUFV1','NDUFV2','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'};
geneNames = {'nuoA','nuoB','nuoC','nuoE','nuoF','nuoG','nuoH','nuoI','nuoJ','nuoK','nuoL','nuoM','nuoN', ...
                                          'appB','appC','atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH','ctaB','ctaC','ctaD','ctaE','ctaF', ...
                                          'cydA','cydB','cydX','cyoA','cyoB','cyoC','cyoD','qcrA', 'qcrB','qcrC','sdhA','sdhB','sdhC','sdhD', ...
'NDUFA9','NDUFS5','NDUFS6','NDUFS4','NDUFB11','NDUFA1','NDUFB10','NDUFA8','NDUFAB1','NDUFB9','NDUFB7','DAP13','GRIM19','NDUFA11','NDUFA6','NDUFA5','NDUFB3','NDUFA2','NDUFA10','NDUFV3','NDUFB2','NDUFB8','NDUFC1','NDUFA4','NDUFB1','NDUFB5','NDUFB6','NDUFB4','NDUFA7','NDUFC2','NDUFA3'};

useTempGenomes = 0;
useHMPGenomes = 1;
count=0;
if useHMPGenomes
   outputDir1 = [outputDir filesep 'runBLASTHMPRefWithSuper'];
   if ~exist(outputDir1,'dir')
       mkdir(outputDir1);
   end
end

parfor z=1:length(geneNames)
    outputDir2 = [outputDir1 filesep geneNames{z}];
    if ~exist(outputDir2,'dir')
        mkdir(outputDir2)
    end
    for j=1:2
        if j==1
            upGenera = ZhangZhaoUpGenera;
            downGenera = ZhangZhaoDownGenera;
            outputDir3 = [outputDir2 filesep 'ZhangZhao'];
        else
            upGenera = ForslundHildebrandUpGenera;
            downGenera = ForslundHildebrandDownGenera;
            outputDir3 = [outputDir2 filesep 'ForslundHildebrand'];
        end
        
        if ~exist(outputDir3,'dir')
            mkdir(outputDir3)
        end
        
        scores = [];
        for k=1:2
            if k==1
                generaToSearch = upGenera;
                outputDir4 = [outputDir3 filesep 'up'];
            else
                generaToSearch = downGenera;
                outputDir4 = [outputDir3 filesep 'down'];
            end
            
            if ~exist(outputDir4,'dir')
                mkdir(outputDir4)
            end

            for i=1:length(generaToSearch)
                %count = count+1;
                if useTempGenomes
                    database = 'nr';
                elseif useHMPGenomes
                    database = ['"/mnt/extra/blast/HMPRef' filesep generaToSearch{i} '"'];
                end            
                commandPrefix = ['blastp -db ' database ' -query /home/ubuntu/MATLAB/GutMicrobiota/input' filesep geneNames{z} '.fasta'];
                if useTempGenomes
                    commandPrefix = [commandPrefix ' -gilist "' outputDir filesep 'extractGIs' filesep generaToSearch{i} '.txt"'];
                end
                command = [commandPrefix ' -outfmt "6 evalue qacc sacc qstart qend sstart send" -out "' outputDir4 filesep generaToSearch{i} '.blast"'];
                disp(command);
                status=system(command);
                %if count > 100
                    %blurt = blurt+1;
                    %end
            end
            
            for i=1:length(generaToSearch)
                FI = fopen([outputDir4 filesep generaToSearch{i} '.blast']);
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