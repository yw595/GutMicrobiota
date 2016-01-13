function runCDDOnce(database,genomeFile,temp,out)

status=system(['rpsblast -db ' database ' -query ' genomeFile ' -evalue 1 -outfmt "6 evalue qacc sacc qstart qend sstart send" -out ' temp]);
status=system(sprintf('sort -k1g < %s > %s',temp,out));
status=system(sprintf('rm %s',temp));

end