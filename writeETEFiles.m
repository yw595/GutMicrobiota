ZhangZhaoGenera = [ZhangZhaoUpGenera; ZhangZhaoDownGenera];
ForslundHildebrandGenera = [ForslundHildebrandUpGenera; ForslundHildebrandDownGenera];

generaTemp = ZhangZhaoGenera;
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
		ZhangZhaoGenera = generaTemp;

outputDir1 = '/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles';
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

FI = fopen('/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles/ZhangZhaoGenera.txt','w');
for i=1:length(ZhangZhaoGenera)
    fprintf(FI,'%s\n',ZhangZhaoGenera{i});
end
fclose(FI);

FI = fopen('/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles/ForslundHildebrandGenera.txt','w');
for i=1:length(ForslundHildebrandGenera)
    fprintf(FI,'%s\n',ForslundHildebrandGenera{i});
end
fclose(FI);
