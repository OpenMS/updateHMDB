from lxml import etree
import sys
import glob
from collections import defaultdict, OrderedDict

missing_value = ""


class OrderedDefaultSetDict(OrderedDict):
    def __missing__(self, key):
        self[key] = value = set()
        return value


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: generateHMDBFilesForOpenMS.py '/path/to/folder/withHMDBxmlfiles/' ")
        exit(1)
    folder = sys.argv[1]
    print("Globbing: " + folder + '*_metabolites.xml')
    files = glob.glob(folder + '*_metabolites.xml')
    print("Reading files: " + str(files))

    if 'hmdb_metabolites.xml' in files and len(files) > 1:
        sys.stderr.write(
            "When running on an old folder layout from HMDB with multiple files, please delete hmdb_metabolites.xml.")
        exit(1)

    version = None
    db2structmapping = open('HMDB2StructMapping.tsv', 'w')
    dbmapping = open('HMDBMappingFile.tsv', 'w')
    formula2id = defaultdict(list)
    #formula2weight = defaultdict(float)
    weight2formula = OrderedDefaultSetDict()
    for x in range(0, len(files)):
        # trigger at end of each important XML element
        context = etree.iterparse(files[x], events=('end',))

        smiles = ""
        inchi = ""
        chemical_formula = ""
        monisotopic_molecular_weight = ""
        accession = None
        name = None

        for event, elem in context:
            # TODO can we set the default namespace in the beginning?
            if version is None and elem.tag == '{http://www.hmdb.ca}version':
                version = elem.text
                continue
            if elem.tag == '{http://www.hmdb.ca}inchi':
                inchi = elem.text
                continue
            if elem.tag == '{http://www.hmdb.ca}chemical_formula':
                chemical_formula = elem.text
                continue
            if elem.tag == '{http://www.hmdb.ca}monisotopic_molecular_weight' and elem.text is not None:
                monisotopic_molecular_weight = elem.text
                continue
            if elem.tag == '{http://www.hmdb.ca}smiles':
                smiles = elem.text
                continue
            if name is None and elem.tag == '{http://www.hmdb.ca}name':
                if elem.text is None:
                    raise RuntimeError("Text of name element returned None")
                name = elem.text.strip()
                continue
            if accession is None and elem.tag == '{http://www.hmdb.ca}accession':
                if elem.text is None:
                    raise RuntimeError("Text of accession element returned None")
                accession = elem.text
                continue
            if elem.tag == '{http://www.hmdb.ca}metabolite':
                if (inchi or smiles) and monisotopic_molecular_weight is not None and chemical_formula != "":
                    # at least one not None or empty string and monoiso weight is available
                    # name and accession can never be None otherwise it fails during parsing of these elements

                    hmdbid = "HMDB:" + accession
                    # db2struct file can be directly written
                    db2structmapping.write(hmdbid + "\t")
                    db2structmapping.write(name + "\t")
                    db2structmapping.write(smiles if smiles is not None else missing_value + "\t")
                    db2structmapping.write(inchi if inchi is not None else missing_value + "\t")
                    db2structmapping.write("\n")

                    # gather mappings
                    #formula2weight[chemical_formula] = float(monisotopic_molecular_weight)
                    weight2formula[float(monisotopic_molecular_weight)].add(chemical_formula)
                    formula2id[chemical_formula].append(hmdbid)
                else:
                    print("Skipping incomplete entry: " + accession)

                # reset variable for next element
                smiles = ""
                inchi = ""
                accession = None
                name = None
                chemical_formula = ""
                monisotopic_molecular_weight = 0.

                # remove the processed element from the tree to save memory (https://lxml.de/3.0/parsing.html)
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]

    # EXTRA metabolites that are added for backwards compatibility
    # https://github.com/OpenMS/OpenMS/pull/1185
    # They seem to have been added to the OpenMS version of the HMDB at some point
    formula2id['C10(2)H3(1)H16NO4'].append("EXTRA:EXTRA001")
    formula2id['C16H18N4O2'].append("EXTRA:EXTRA002")
    formula2id['C12(2)H6(1)H8N4O4S'].append("EXTRA:EXTRA003")
    formula2id['C33H40N2O9'].append("EXTRA:EXTRA004")
    formula2id['C32H41NO2'].append("EXTRA:EXTRA005")
    formula2id['C23(2)H3(1)H42NO4'].append("EXTRA:EXTRA006")
    formula2id['C25(2)H3(1)H46NO4'].append("EXTRA:EXTRA007")

    #formula2weight['C10(2)H3(1)H16NO4'] = 220.1502383412
    #formula2weight['C16H18N4O2'] = 298.1429758428
    #formula2weight['C12(2)H6(1)H8N4O4S'] = 316.111236124
    #formula2weight['C23(2)H3(1)H42NO4'] = 402.3536891758
    #formula2weight['C25(2)H3(1)H46NO4'] = 430.3849893042

    weight2formula[220.1502383412].add('C10(2)H3(1)H16NO4')
    weight2formula[298.1429758428].add('C16H18N4O2')
    weight2formula[316.111236124].add('C12(2)H6(1)H8N4O4S')
    weight2formula[402.3536891758].add('C23(2)H3(1)H42NO4')
    weight2formula[430.3849893042].add('C25(2)H3(1)H46NO4')

    #formula2mass_sortedbymass = sorted(formula2weight.items(), key=operator.itemgetter(1))

    # Write header of HMDBMappingFile.tsv
    dbmapping.write('database_name')
    dbmapping.write('\t')
    dbmapping.write('HMDB')
    dbmapping.write('\n')
    dbmapping.write('database_version')
    dbmapping.write('\t')
    dbmapping.write(version)
    dbmapping.write('\n')

    # Write rows:
    for weight, formulas in weight2formula.items():
        if len(formulas) > 1:
            print("WARNING: two equal formulas with different masses found")
        for formula in formulas:
            if not formulas:
                continue
            dbmapping.write(str(weight))
            dbmapping.write('\t')
            dbmapping.write(formula)

            for x in formula2id[formula]:
                dbmapping.write('\t')
                dbmapping.write(x)
            dbmapping.write('\n')


if __name__ == '__main__':
    main()
