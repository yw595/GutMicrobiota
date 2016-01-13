function [score matchName qCovS pIdent blastOcc] = readBLASTFile(BLASTFile,dispLines)

    if ~exist('dispLines','var')
        dispLines=1;
    end
    
FI = fopen(BLASTFile);
blastOcc=0;
if FI ~= -1
    line = fgetl(FI);
    if line ~= -1
        words = strsplit(line,'\t');
        regExpString = '(putative|hypothetical|uncharacterized|predicted)';
        while line ~= -1 & ~isempty(regexp(words{9}, regExpString))
            words = strsplit(line,'\t');
            line = fgetl(FI);
            if dispLines
                disp(line)
            end
        end
        if isempty(regexp(words{9},regExpString))
            score = str2num(words{1});
            matchName = words{9};
            qCovS = str2num(words{10});
            pIdent = str2num(words{11});
            blastOcc = score<.00001 && qCovS>.3 && pIdent>.3;
        else
            score = -1;
            matchName = 'NA';
            qCovS = -1;
            pIdent = -1;
        end
    else
        score = -1;
        matchName = 'NA';
        qCovS = -1;
        pIdent = -1;
    end
    fclose(FI);
else
    score = -2;
    matchName = 'NA';
    qCovS = -2;
    pIdent = -2;
end

end