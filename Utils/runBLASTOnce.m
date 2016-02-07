function runBLASTOnce(database,query,gilist,out)
commandPrefix = ['blastp -db ' database ' -query ' query];
if ~strcmp(gilist,'')
    commandPrefix = [commandPrefix ' -gilist ' gilist];
end
command = [commandPrefix [' -outfmt "6 evalue qacc sacc qstart qend sstart send qseqid sseqid qcovs pident" -out '] out];
disp(command);
status = system(command);
end