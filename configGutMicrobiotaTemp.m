outputDir = '/home/ubuntu/MATLAB/GutMicrobiota/output';
HMPDir = '/mnt/extra/blast/HMPRef';
tempDir = '/mnt/extra/blast/CDDGenomes/temp';
CDDDir = '/mnt/extra/blast/CDD';
inputDir = '/home/ubuntu/MATLAB/GutMicrobiota/input';

complexesToGenes = containers.Map({'mitochondrial_complex_I_supernumerary','mitochondrial_complex_I','ndh','Complex I','cytochrome bd-II oxidase','ATPase','B. sub cytochrome c oxidase','cytochrome bd-I oxidase','cytochrome bo oxidase','B. sub cytochrome b/c1 oxidase','Complex II','butyrate_fermentation','folate_cycle'}, ...
                                      {{'NDUFA9','NDUFS5','NDUFS6','NDUFS4','NDUFB11','NDUFA1','NDUFB10','NDUFA8','NDUFAB1','NDUFB9','NDUFB7','DAP13','GRIM19','NDUFA11','NDUFA6','NDUFA5','NDUFB3','NDUFA2','NDUFA10','NDUFV3','NDUFB2','NDUFB8','NDUFC1','NDUFA4','NDUFB1','NDUFB5','NDUFB6','NDUFB4','NDUFA7','NDUFC2','NDUFA3'}, ...
                        {'NDUFS1','NDUFS2','NDUFS3','NDUFS7','NDUFS8', 'NDUFV1','NDUFV2','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'}, ...
                   {'ndh'},{'nuoA','nuoB','nuoC','nuoE','nuoF','nuoG','nuoH','nuoI','nuoJ','nuoK','nuoL','nuoM','nuoN'}, ...
    {'appB','appC'},{'atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH'},{'ctaC','ctaD','ctaE','ctaF'}, ...
    {'cydA','cydB'}, {'cyoA','cyoB','cyoC','cyoD'}, {'qcrA','qcrB','qcrC'}, {'sdhA','sdhB','sdhC','sdhD'}, ...
    {'thiA','hbd','crt','bcd','ptb','rnfA','rnfB','rnfC','rnfD','rnfE'},{'MS','MAT','AHCY','MTHFR','TS','DHFR','FTHFDH','MTCH_MTHD','SHMT','GDC','CGS','CBL'}});

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
            words = strsplit(line,'\t');
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
            complexCdds{end+1} = cddName;
            line = fgetl(FI);
        end
        complexesToCdds(complexName) = complexCdds;
        fclose(FI);
    end
end

complexesToGenesAnnot=containers.Map('NADH_dehydrogenase', 'cytochrome_bd_oxidase','cytochrome_bo_oxidase','cytochrome_c_oxidase','cytochrome_d_oxidase','cytochrome_o_oxidase','ATP_synthase','succinate_dehydrogenase'},{{'A','B','C','D','E','F','G','H','I','J','K','L','M','N'},{'I_I','I_II','II_I','II_II'},{''},{'I','II','III','IV'},{'I','II'},{'I','II','III','IV'},{'alpha','beta','gamma','delta','epsilon','A','B','C'},{''}});

complexesToComplexes = containers.Map({'Complex I','ATPase','B. sub cytochrome c oxidase','B. sub cytochrome b/c1 oxidase','Complex II'},{'NADH_dehydrogenase','ATP_synthase','cytochrome_c_oxidase','cytochrome_c_oxidase','succinate_dehydrogenase'});
complexes = keys(complexesToComplexes);

genesToGenes = containers.Map();
for i=1:length(complexes)
    if strcmp(complexes{i},'Complex I')
        localGeneMap = containers.Map({'G',{'C','D'},{'C','D'},'B','I','F','E','H','N','A','M','K','L','J'},{'NDUFS1','NDUFS2','NDUFS3','NDUFS7','NDUFS8', 'NDUFV1','NDUFV2','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'});
    elseif strcmp(complexes{i},'cytochrome bd-II oxidase')
        localGeneMap = containers.Map({'appC','appB'},{'II_I','II_II'});
    elseif strcmp(complexes{i},'cytochrome bd-I oxidase')
        localGeneMap = containers.Map({'cydA','cydB'},{'I_I','I_II'});
    elseif strcmp(complexes{i},'cytochrome bo oxidase')
        localGeneMap = containers.Map({'cyoA','cyoB','cyoC','cyoD'},{'II','I','III','IV'});
    elseif strcmp(complexes{i},'B. sub cytochrome c oxidase')
        localGeneMap = containers.Map({'ctaD','ctaC','ctaE','ctaF'},{'I','II','III','IV'});
    elseif strcmp(complexes{i},'ATPase')
        localGeneMap = containers.Map({'alpha','A','epsilon','beta','C','B','gamma','delta'},{'atpA','atpB','atpC','atpD','atpE','atpF','atpG','atpH'});
        elseif strcmp()
    end
end

readZhangZhao;
readForslundHildebrand;