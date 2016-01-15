from ete3 import NCBITaxa
import itertools

ZhangZhaoGenera = []
with open('/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles/ZhangZhaoGenera.txt') as f:
    for line in f:
        ZhangZhaoGenera.append(line.strip())
    f.close()

ForslundHildebrandGenera = []
with open('/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles/ForslundHildebrandGenera.txt') as f:
    for line in f:
        ForslundHildebrandGenera.append(line.strip())
    f.close()

RefSeqSpecies = []
with open('/home/ubuntu/MATLAB/GutMicrobiota/input/reference_genomes.txt') as f:
    next(f)
    for line in f:
        words = line.split('\t')
        RefSeqSpecies.append(words[0])
    f.close()

ncbi = NCBITaxa()
name2taxid = ncbi.get_name_translator(list(set(ZhangZhaoGenera+ForslundHildebrandGenera+RefSeqSpecies)))

tree = ncbi.get_topology(list(itertools.chain.from_iterable(list(name2taxid.values()))),intermediate_nodes=True)
#print(tree.get_ascii(attributes=['sci_name']), file=open('/home/ubuntu/taxonomy.txt','w'))

#print(tree.name)
# fh = open('/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles/closestSpecies.txt','w')
# for genus in ZhangZhaoGenera+ForslundHildebrandGenera:
#     print(genus)
#     minDist = -1
#     minDistSpecies = ''
#     for species in RefSeqSpecies:
#         genusNode = tree.search_nodes(name=str(name2taxid[genus][0]))[0]
#         speciesNode = tree.search_nodes(name=str(name2taxid[species][0]))[0]
#         dist = tree.get_distance(speciesNode, genusNode)
#         if minDist == -1:
#             minDist = dist
#             minDistSpecies = species
#         elif minDist > dist:
#             minDist = dist
#             minDistSpecies = species
#     print(genus+" "+minDistSpecies+" "+str(minDist))
#     print(genus+"\t"+minDistSpecies,file=fh)
# fh.close()

fh = open('/home/ubuntu/MATLAB/GutMicrobiota/output/writeETEFiles/allDescendants.txt','w')
for genus in ZhangZhaoGenera+ForslundHildebrandGenera:
    print(genus)
    genusNode = tree.search_nodes(name=str(name2taxid[genus][0]))[0]
    #descendants = genusNode.get_descendants()
    #descendantNames = []
    #for d in descendants:
    #    descendantNames.append(d.name)
    descendants = ncbi.get_descendant_taxa(str(name2taxid[genus][0]),intermediate_nodes=True)
    descendantNames = [str(name2taxid[genus][0])]
    for d in descendants:
        descendantNames.append(str(d))
    print(genus+"\t"+",".join(descendantNames),file=fh)
fh.close()
