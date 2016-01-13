function matrix = indexPar(inMatrix,val,indices)

matrix=inMatrix;
matrix(indices(1),indices(2),indices(3),indices(4),indices(5))=val;

end