configGutMicrobiota;
inputDir = '/home/ubuntu/MATLAB/GutMicrobiota/output/extractGIs';
outputDir1 = [outputDir filesep 'runCDD'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end


ZhangZhaoUpGenera = changeNames(ZhangZhaoUpGenera);
ZhangZhaoDownGenera = changeNames(ZhangZhaoDownGenera);
ForslundhildebrandUpGenera = changeNames(ForslundHildebrandUpGenera);
ForslundHildebrandDownGenera = changeNames(ForslundHildebrandDownGenera);

files = dir(inputDir);
base1 = 'wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=';
base2 = '&rettype=fasta"';
for i=1:length(files)
    if ~isempty(regexp(files(i).name,'.txt'))
        %species = files(i).name(1:regexp(files(i).name,'.txt')-1);
        %disp(species);
        
        % FI = fopen([inputDir filesep files(i).name]);
        % dataFields = textscan(FI,'%s');
        % dataFields = dataFields{1};
        % count = 0;
        % while count*400+1 <= length(dataFields)
        %     disp(count)
        %     payload = strjoin(dataFields(count*400+1:min(length(dataFields),count*400+400)),',');
        %     if count==0
        %         [status cmdout]=system([base1 payload base2 ' -O - >' outputDir1 filesep species '.temp.fasta']);
        %     else                
        %         [status cmdout]=system([base1 payload base2 ' -O - >>' outputDir1 filesep species '.temp.fasta']);
        %     end
        %     if isempty(regexp(cmdout,'Name or service not known'))
        %         count = count+1;
        %     end
        % end
        % status = system(['grep -v "^$" ' outputDir1 filesep species '.temp.fasta > ' outputDir1 filesep species '.fasta']);
        % status = system(['rm ' outputDir1 filesep species '.temp.fasta']);
        % fclose(FI);
    end
end

useTempGenomes = 0;
useHMPGenomes = 1;
if useHMPGenomes
   outputDir1 = [outputDir filesep 'runCDDHMPRef'];
   if ~exist(outputDir1,'dir')
       mkdir(outputDir1);
   end
end
if useTempGenomes
    CDDGenomesDir = '/mnt/extra/blast/CDDGenomes/temp';
elseif useHMPGenomes
    CDDGenomesDir = '/mnt/extra/blast/HMPRef';
end
%fastaFiles = dir(CDDGenomesDir);
CDDDir = '/mnt/extra/blast/CDD';
complexFiles = {'NADH_dehydrogenase_CDD.temp.txt','cytochrome_CDD.temp.txt','ATP_synthase_CDD.temp.txt'};
count = 0;
cddIDs = {};
cddNames = {};
for i=1:length(complexFiles)
    FI1 = fopen([CDDDir filesep complexFiles{i}]);
    line = fgetl(FI1);
    while line~=-1
        words = strsplit(line,'\t');
        cddGi = words{1};
        cddID = words{2};
        cddIDs{end+1} = cddID;
        cddName = words{4};
        cddName = strrep(cddName,' ','_');
        cddName = strrep(cddName,';','');
        cddName = strrep(cddName,'\.','');
        cddName = [cddID '_' cddName];
        cddNames{end+1} = cddName;
        line = fgetl(FI1);
        
        % FI2 = fopen(['/mnt/extra/blast/CDD/smp' filesep cddID '.list'],'w');
        % fprintf(FI2,'%s\n',[cddID '.smp']);
        % fclose(FI2);
        % cd('/mnt/extra/blast/CDD/smp');
        % status=system(['makeprofiledb -in ' cddID '.list']);
        % cd('/home/ubuntu/MATLAB')
    end
    fclose(FI1);
end

parfor z=1:length(cddNames)
    cddName = cddNames{z}; cddID = cddIDs{z};
    disp(cddID)
    outputDir2 = [outputDir1 filesep cddName];
    if ~exist(outputDir2,'dir')
        mkdir(outputDir2);
    end

    %FI2 = fopen([outputDir1 filesep 'temp.txt'],'w');
    %fprintf(FI2,'%s\n',cddGi);
    %fclose(FI2);
    for j=1:2
        if j==1
            upGenera = ForslundHildebrandUpGenera;
            downGenera = ForslundHildebrandDownGenera;
            outputDir3 = [outputDir2 filesep 'ForslundHildebrand'];
        else
            upGenera = ZhangZhaoUpGenera;
            downGenera = ZhangZhaoDownGenera;
            outputDir3 = [outputDir2 filesep 'ZhangZhao'];
        end
        if ~exist(outputDir3,'dir')
            mkdir(outputDir3);
        end

        for k=1:2
            if k==1
                genera = upGenera;
                outputDir4 = [outputDir3 filesep 'up'];
            else
                genera = downGenera;
                outputDir4 = [outputDir3 filesep 'down'];
            end
            if ~exist(outputDir4,'dir')
                mkdir(outputDir4);
            end             

            for l=1:length(genera)
                genera{l} = strrep(genera{l},' ','_');
                %count = count+1;
                %count
                if useTempGenomes
                    genomeFile = [CDDGenomesDir filesep genera{l} '.test.fasta'];
                elseif useHMPGenomes
                    genomeFile = [CDDGenomesDir filesep genera{l} '.faa'];
                    if ~exist(genomeFile,'file')
                        genomeFile = ['/mnt/extra/blast/CDDGenomes/temp' filesep genera{l} '.test.fasta'];
                    end
                end
                status=system(['rpsblast -db ' CDDDir  filesep 'smp' filesep cddID ' -query ' genomeFile ' -evalue 1 -outfmt "6 evalue qacc sacc qstart qend sstart send" -out ' outputDir4 filesep 'temp.blast']);
                status=system(sprintf('sort -k1g < %s > %s',[outputDir4 filesep 'temp.blast'],[outputDir4 filesep genera{l} '.blast']));
                status=system(sprintf('rm %s',[outputDir4 filesep 'temp.blast']));
                %if count > 10
                %    blurt = blurt+1;
                %end                    
            end
        end
    end
    %status=system(['rm ' outputDir1 filesep 'temp.txt']);
end

% parfor i=1:length(cddIDs)
%     cddID = cddIDs{i};
%     FI2 = fopen(['/mnt/extra/blast/CDD/smp' filesep cddID '.list'],'w');
%     fprintf(FI2,'%s\n',[cddID '.smp']);
%     fclose(FI2);
%     cd('/mnt/extra/blast/CDD/smp');
%   o  status=system(['makeprofiledb -in ' cddID '.list -out ' cddID ' -title ' cddID]);
%     cd('/home/ubuntu/MATLAB')
% end
