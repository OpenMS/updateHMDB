from lxml import etree
import os
import csv
from os import getcwd
from os import listdir
import sys
import re

files = os.listdir('hmdb_metabolites')
#remove .DS_Store
if (files[0]=='.DS_Store'):
    files.pop(0)

if 'hmdb_metabolites.xml' in files and len(files)>1:
    print "please"
    sys.exit()

f = open('HMDB2StructMapping.tsv', 'w')
for x in range(0,len(files)):
#for x in range(0,10):
   context = etree.iterparse('hmdb_metabolites/'+ files[x], events=('end',))

   smiles = ""
   inchi = ""
   accession = None

   for action, elem in context:
      if elem.tag=='{http://www.hmdb.ca}inchi':
	 inchi = elem.text
         pass
      if elem.tag=='{http://www.hmdb.ca}smiles':
	 smiles = elem.text
         pass
      if elem.tag=='{http://www.hmdb.ca}name':
	 name = elem.text.strip()
         pass
      if accession==None and elem.tag=='{http://www.hmdb.ca}accession':
	 accession = elem.text
         pass
      if elem.tag=='{http://www.hmdb.ca}metabolite':
         if(inchi or smiles):
            f.write("HMDB:"+accession + "\t")
            f.write(name.encode('utf-8') + "\t")
            f.write(smiles + "\t")
            f.write(inchi + "\t")
            f.write("\n")
         smiles = ""
         inchi = ""
         accession = None
         elem.clear()
         while elem.getprevious() is not None:
            del elem.getparent()[0]


f2 = open('extras.tsv', 'r')
extras = f2.readlines()
for line in extras:
	elements = re.split(r'\t+', line)
	for element in elements[:-1]:
		f.write(element)
		f.write('\t')
	f.write(elements[-1])

f.close()
