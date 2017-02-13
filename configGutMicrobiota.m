extraDir = '/mnt/extra/blast/HMPRef';
fileNames = dir(extraDir);
generaToFileNames = containers.Map;
for i=1:length(fileNames)
    % INSIDIOUS BUG, without underscore, will match genusName.faa
    % from previous runs, then genusName variable will be empty,
    % resulting in bunch of fileNames linked to ''
    if ~isempty(regexp(fileNames(i).name,'_.*\.faa$'))
        underscoreIdx = regexp(fileNames(i).name,'_');
        genusName = fileNames(i).name(1:underscoreIdx-1);
        %disp(genusName)
        if isKey(generaToFileNames,genusName)
            genusFileNames = generaToFileNames(genusName);
            genusFileNames{end+1} = fileNames(i).name;
            generaToFileNames(genusName) = genusFileNames;
        else
            generaToFileNames(genusName) = {fileNames(i).name};
        end
    end
end

outputDir = '/home/ubuntu/MATLAB/GutMicrobiota/output';
HMPDir = '/mnt/extra/blast/HMPRef';
tempDir = '/mnt/extra/blast/CDDGenomes/temp';
CDDDir = '/mnt/extra/blast/CDD';
inputDir = '/home/ubuntu/MATLAB/GutMicrobiota/input';

complexes = {'Complex_I','Complex_II','B_sub_cytochrome_bc1_oxidase','B_sub_cytochrome_c_oxidase','cytochrome_bo_oxidase','cytochrome_bd_II_oxidase','cytochrome_bd_I_oxidase','ATP_synthase','ndh','mitochondrial_complex_I','mitochondrial_complex_I_supernumerary','butyrate_fermentation','folate_cycle'};
%complexesToGenes =containers.Map({'mitochondrial_complex_I_supernumerary','mitochondrial_complex_I','ndh','Complex_I','cytochrome_bd_II_oxidase','ATP_synthase','B_sub_cytochrome_c_oxidase','cytochrome_bd_I_oxidase','cytochrome_bo_oxidase','B_sub_cytochrome_bc1_oxidase','Complex_II','butyrate_fermentation','folate_cycle'},...
complexesToGenes = containers.Map(complexes, ...
                                      {{'nuoA','nuoB','nuoC','nuoE','nuoF','nuoG','nuoH','nuoI','nuoJ','nuoK','nuoL','nuoM','nuoN'},{'sdhA','sdhB','sdhC','sdhD'},{'qcrA','qcrB','qcrC'},{'ctaC','ctaD','ctaE','ctaF'},{'cyoA','cyoB','cyoC','cyoD'},{'appB','appC'},{'cydA','cydB'},{'atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH'},{'ndh'},{'NDUFS1','NDUFS2','NDUFS3','NDUFS7','NDUFS8', 'NDUFV1','NDUFV2','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'},{'NDUFA9','NDUFS5','NDUFS6','NDUFS4','NDUFB11','NDUFA1','NDUFB10','NDUFA8','NDUFAB1','NDUFB9','NDUFB7','DAP13','GRIM19','NDUFA11','NDUFA6','NDUFA5','NDUFB3','NDUFA2','NDUFA10','NDUFV3','NDUFB2','NDUFB8','NDUFC1','NDUFA4','NDUFB1','NDUFB5','NDUFB6','NDUFB4','NDUFA7','NDUFC2','NDUFA3'},{'thiA','hbd','crt','bcd','ptb','rnfA','rnfB','rnfC','rnfD','rnfE'},{'MS','MAT','AHCY','MTHFR','TS','DHFR','FTHFDH','MTCH_MTHD','SHMT','GDC','CGS','CBL'}});

if ~exist('cddIDs','var')
    complexFiles = {'NADH_dehydrogenase_CDD.temp.txt','cytochrome_CDD.temp.txt','ATP_synthase_CDD.temp.txt'};
    cddIDs = {};
    cddNames = {};
    complexesToCdds = containers.Map;
    for i=1:length(complexFiles)
        FI = fopen([CDDDir filesep complexFiles{i}]);
        complexName = complexFiles{i};
        complexName = complexName(1:regexp(complexName,'_CDD')-1);
        complexCdds = {};
        line = fgetl(FI);
        while line~=-1
            words = strsplitYiping(line,'\t');
            cddGi = words{1};
            cddID = words{2};
            cddIDs{end+1} = cddID;
            cddName = words{4};
            cddName = strrep(cddName,' ','_');
            cddName = strrep(cddName,';','');
            cddName = strrep(cddName,'\.','');
            cddName = strrep(cddName,'''','_prime_');
            cddName = [cddID '_' cddName];
            cddNames{end+1} = cddName;
            complexCdds{end+1} = [cddName ' ' cddID];
            line = fgetl(FI);
        end
        complexesToCdds(complexName) = complexCdds;
        fclose(FI);
    end
end

%complexesToGenesAnnot=containers.Map({'NADH_dehydrogenase', 'cytochrome_bd_oxidase','cytochrome_bo_oxidase','cytochrome_c_oxidase','ATP_synthase','succinate_dehydrogenase','cytochrome_bc1_oxidase'},{{'A','B','C','D','E','F','G','H','I','J','K','L','M','N'},{'I_I','I_II','II_I','II_II'},{'I','II','III','IV'},{'I','II','III','IV','VI','cbb3_I','cbb3_II','cbb3_III','cbb3_IV'},{'alpha','beta','gamma','delta','epsilon','A','B','C'},{'A','B','C','D'},{'B','C1','ISP'}});
complexesAnnot = {'NADH_dehydrogenase','succinate_dehydrogenase','cytochrome_bc1_oxidase','cytochrome_c_oxidase','cytochrome_bo_oxidase','cytochrome_bd_oxidase','ATP_synthase'};
complexesToGenesAnnot=containers.Map(complexesAnnot,{{'A','B','C','D','E','F','G','H','I','J','K','L','M','N'},{'A','B','C','D'},{'B','C1','ISP'},{'I','II','III','IV'},{'I','II','III','IV'},{'I_I','I_II','II_I','II_II'},{'alpha','beta','gamma','delta','epsilon','A','B','C'}});

matchingComplexes = {'Complex_I','Complex_II','B_sub_cytochrome_bc1_oxidase','B_sub_cytochrome_c_oxidase','cytochrome_bo_oxidase','cytochrome_bd_II_oxidase','cytochrome_bd_I_oxidase','ATP_synthase'};
complexesToComplexes = containers.Map(matchingComplexes,{'NADH_dehydrogenase','succinate_dehydrogenase','cytochrome_bc1_oxidase','cytochrome_c_oxidase','cytochrome_bo_oxidase','cytochrome_bd_oxidase','cytochrome_bd_oxidase','ATP_synthase'});

genesToGenes = containers.Map();
for i=1:length(matchingComplexes)
    if strcmp(matchingComplexes{i},'Complex_I')
        localGeneMap = containers.Map({'G','C','D','B','I','F','E','H','N','A','M','K','L','J'},{'NDUFS1','NDUFS2','NDUFS3','NDUFS7','NDUFS8', 'NDUFV1','NDUFV2','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'});
    elseif strcmp(matchingComplexes{i},'cytochrome_bd_II_oxidase')
        localGeneMap = containers.Map({'II_I','II_II'},{'appC','appB'});
    elseif strcmp(matchingComplexes{i},'cytochrome_bd_I_oxidase')
        localGeneMap = containers.Map({'I_I','I_II'},{'cydA','cydB'});
    elseif strcmp(matchingComplexes{i},'cytochrome_bo_oxidase')
        localGeneMap = containers.Map({'II','I','III','IV'},{'cyoA','cyoB','cyoC','cyoD'});
    elseif strcmp(matchingComplexes{i},'B_sub_cytochrome_c_oxidase')
        localGeneMap = containers.Map({'I','II','III','IV'},{'ctaD','ctaC','ctaE','ctaF'});
    elseif strcmp(matchingComplexes{i},'ATP_synthase')
        localGeneMap = containers.Map({'alpha','A','epsilon','beta','C','B','gamma','delta'},{'atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH'});
    elseif strcmp(matchingComplexes{i},'Complex_II')
        localGeneMap = containers.Map({'A','B','C','D'},{'sdhA','sdhB','sdhC','sdhD'});
    elseif strcmp(matchingComplexes{i},'B_sub_cytochrome_bc1_oxidase')
        localGeneMap = containers.Map({'ISP','B','C1'},{'qcrA','qcrB','qcrC'});
    else
        disp('ERROR: No localGeneMap');
        disp(matchingComplexes{i});
        break;
    end

    localGenes = keys(localGeneMap);
    for j=1:length(localGenes)
        genesToGenes([complexesToComplexes(matchingComplexes{i}) '_' localGenes{j}]) = localGeneMap(localGenes{j});
    end
end

matchingGenesArray = {};
for i=1:length(complexesAnnot)
    geneNamesAnnot = complexesToGenesAnnot(complexesAnnot{i});
    for j=1:length(geneNamesAnnot)
        matchingGenesArray{end+1} = genesToGenes([complexesAnnot{i} '_' geneNamesAnnot{j}]);
    end
end


readZhangZhao;
readForslundHildebrand;
uniqGenera = unique([ZhangZhaoUpGenera ZhangZhaoDownGenera ForslundHildebrandUpGenera ForslundHildebrandDownGenera]);
avoidGenera={'Allobaculum','Clostridium','Butyricimonas', 'Oscillibacter','candidate_division_TM7','Rothia','Shigella'};
uniqGeneraAnnot = uniqGenera(~ismember(uniqGenera,avoidGenera));