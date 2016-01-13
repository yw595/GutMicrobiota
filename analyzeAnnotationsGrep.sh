#!/bin/bash

grep -i -P "(ATP.*syn|NADH.*dehydro|ubiquin.*dehydro|succin.*dehydro|cytochrome)" *.faa | grep -i -P -v "(biosynthesis|biogenesis|assembly|maturation|flagellum|P450|V-type|putative|lactate.*dehydro|butan.*dehydro|fumarate.*dehydro|glucon.*dehydro|formate.*dehydro|lipoyl|semialdehyde.*dehydro|thiosulfate)" | cut -d" " -f1 --complement | grep -i -P -v "\.faa" | sort | uniq > uniqSortedAllETC.txt
grep -i -P "(NADH.*reduct|ubiquin.*reduct|succin.*reduct)" *.faa | cut -d" " -f1 --complement | grep -i -P -v "\.faa" | sort | uniq > uniqReduct.txt
grep -P "Rieske" *.faa > rieske.txt
grep -P "menaquin" *.faa | grep -i -P -v "(biosynthesis|biogenesis|assembly|maturation|fumarate|uncharacterized|putative|hypothetical|metabolism|predicted|aryl)" > menaquin.txt
grep -P "Fe-*S" *.faa | grep -i -P -v "(biosynthesis|biogenesis|assembly|maturation|fumarate|uncharacterized|putative|hypothetical|metabolism|predicted|aryl)" > FeS.txt
grep -P -i "(Complex I|Complex V)" *.faa > complex.txt
