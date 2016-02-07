readForslundHildebrand;
readZhangZhao;

%order uncertain
system('python /home/ubuntu/MATLAB/GutMicrobiota/makeETEFiles.py');
writeETEFiles;
extractGIs;

getRefSeq;
readHMPRef;
makeHMPRefBLAST;
system('/home/ubuntu/MATLAB/GutMicrobiota/analyzeAnnotationsGrep.sh');
analyzeAnnotations;

runSeqAnalyses;
makePNGs;