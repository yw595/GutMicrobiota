configGutMicrobiota
outputDir1 = [outputDir filesep 'extractGIs'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

extraDir = '/mnt/extra/blast';
%generaToSearch = [ZhangZhaoUpGenera; ZhangZhaoDownGenera; ForslundHildebrandUpGenera; ForslundHildebrandDownGenera];
namesDmp = [extraDir filesep 'names.dmp'];
namesDmpFind = [outputDir1 filesep 'names.dmp.find'];
% command = 'grep -P "\t(';
% for i=1:length(generaToSearch)-1
%     command = [command generaToSearch{i} '|'];
% end
% command = [command generaToSearch{end} ')\t" ' namesDmp ' | cut -f3,1 > ' namesDmpFind];
% disp(command)
% status=system(command);

FI = fopen([outputDir filesep 'writeETEFiles' filesep 'allDescendants.txt']);
FI1 = fopen(namesDmpFind,'w');
line = fgetl(FI);
while line~=-1
    words = strsplit(line,'\t');
    words1 = strsplit(words{2},',');
    for i=1:length(words1)
        fprintf(FI1,'%s\t%s\n',words1{i},words{1});
    end
    line = fgetl(FI);
end
fclose(FI); fclose(FI1);

FI = fopen(namesDmpFind);
line = fgetl(FI);
namesToTaxids = containers.Map;
while line ~= -1
    words = strsplit(line,'\t');
    namesToTaxids(words{2}) = words{1};
    line = fgetl(FI);
end
fclose(FI);

command = ['gawk ''BEGIN {FS=OFS="\t"}; FNR==NR{array[$1]=$2; next}{if($2 in array){print $1 > "' outputDir1 filesep '" array[$2] ".txt"}}'' ' outputDir1 filesep 'names.dmp.find ' extraDir filesep 'gi_taxid_prot.dmp'];
disp(command);
status=system(command);