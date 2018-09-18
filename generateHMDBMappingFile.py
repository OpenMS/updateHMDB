from xml.dom import minidom
from lxml import etree
import os
import csv
from os import getcwd
from os import listdir
import sys

import operator

files = os.listdir('hmdb_metabolites')
#remove .DS_Store
if (files[0]=='.DS_Store'):
    files.pop(0)

if 'hmdb_metabolites.xml' in files and len(files)>1:
    print "please"
    sys.exit()

f = open('HMDBMappingFile.tsv', 'w')
mapping = dict()
mapping2 = dict()
for x in range(0,len(files)):
#for x in range(0,100):
   context = etree.iterparse('hmdb_metabolites/'+ files[x], events=('end',))

   chemical_formula = "no_entry"
   monisotopic_moleculate_weight = 0
   accession = None

   for action, elem in context:
      if elem.tag=='{http://www.hmdb.ca}chemical_formula':
         chemical_formula = elem.text
         pass
      if elem.tag=='{http://www.hmdb.ca}monisotopic_molecular_weight': # and elem.text!=None:
         monisotopic_moleculate_weight = elem.text
         pass
      if accession==None and elem.tag=='{http://www.hmdb.ca}accession':
         accession = elem.text
         pass
      if elem.tag=='{http://www.hmdb.ca}metabolite':
         #print accession, monisotopic_moleculate_weight
         if monisotopic_moleculate_weight!=None:
            mapping2[chemical_formula] = float(monisotopic_moleculate_weight)
            if not chemical_formula in mapping:
               mapping[chemical_formula] = ["HMDB:" + accession]
            else:
               mapping[chemical_formula].append('HMDB:'+accession)
         chemical_formula = "no_entry"
         monisotopic_moleculate_weight = 0
         accession = None
         elem.clear()
         while elem.getprevious() is not None:
            del elem.getparent()[0]

mapping['C10(2)H3(1)H16NO4'] = ["EXTRA:EXTRA001"]
mapping['C16H18N4O2'] = ["EXTRA:EXTRA002"]
mapping['C12(2)H6(1)H8N4O4S'] = ["EXTRA:EXTRA003"]
mapping['C23(2)H3(1)H42NO4'] = ["EXTRA:EXTRA006"]
mapping['C25(2)H3(1)H46NO4'] = ["EXTRA:EXTRA007"]
mapping2['C10(2)H3(1)H16NO4'] = 220.1502383412
mapping2['C16H18N4O2'] = 298.1429758428
mapping2['C12(2)H6(1)H8N4O4S'] = 316.111236124
mapping2['C23(2)H3(1)H42NO4'] = 402.3536891758
mapping2['C25(2)H3(1)H46NO4'] = 430.3849893042
mapping['C32H41NO2'].append("EXTRA:EXTRA005")
mapping['C33H40N2O9'].append("EXTRA:EXTRA004")
sorted_x = sorted(mapping2.items(), key=operator.itemgetter(1))

f.write('database_name')
f.write('\t')
f.write('HMDB')
f.write('\n')
f.write('database_version')
f.write('\t')
f.write('4.0')
f.write('\n')

for (key,weight) in sorted_x:
   if key=="no_entry":
      continue
   f.write(str(weight))
   f.write('\t')
   f.write(key)
   for x in mapping[key]:
      f.write('\t')
      f.write(x)
   f.write('\n')
