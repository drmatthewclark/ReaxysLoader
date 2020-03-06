# rmc file release number
dbname='mclark'
debug = False

import xml.etree.ElementTree as ET
import psycopg2 as psql
from   psycopg2 import sql
from psycopg2.extensions import AsIs
import glob
import gzip

# for rdkit Mol to Smiles
from rdkit import Chem
from rdkit import RDLogger
import rdkit.Chem.rdChemReactions



def readnextRDfile(file):
    line = '' # file.readline()
    rdfile = ''
    tags = {}
    blankcount = 0

    line = file.readline().rstrip()
    if line == '':
        return None

    while not line.startswith('$RFMT $RIREG'):
        line = file.readline()
        if not line.startswith('$RFMT $RIREG'):
            rdfile = rdfile + line
        # python file reading is brain dead and returns an empty string
        # at the EOF instad of  None or any other EOF marker so we have to
        # guess that a bunch of blank lines means that we have reached the
        # end of file 
        if line == '':
            blankcount += 1
        if blankcount > 2:
            break

    tags = processRXN(rdfile)
    return tags # dictionary 


global counter 
counter = 0

def processRXN(rdfile):
    global counter
    header = '' + '\n' + 'GSMACCS-II07189510252D 1   0.00366     0.00000     0' + '\n\n'  + '  0  0  0     0  0            999 V3000' + '\n'
    tail = 'M  END' 
    counter += 1
    reactants = list()
    products = list()
    currentstructure = header
    state = ''
    dtype = None
    data = {}

    for line in rdfile.splitlines():
        endctab = None
        if line.startswith('M  V30 BEGIN REACTANT'): # reactant section
            state = 'reactant'  # the structures are reactants
            continue
        if line.startswith('M  V30 END REACTANT'):
            state = None   # reactant section has ended
            continue
        if line.startswith('M  V30 BEGIN PRODUCT'):
            state = 'product'  # molecules are now products
            continue
        if line.startswith('M  V30 END PRODUCT'):
            state = None  # product section has ended
            continue
        if line.startswith('M  V30 END CTAB'):
            endctab = True 
        if state:
            currentstructure += line + '\n'
       
        if endctab:
            endctab = None
            if state == 'reactant':
                currentstructure += tail
                reactants.append(currentstructure)
                currentstructure = header
            elif state == 'product':
                currentstructure += tail
                products.append(currentstructure) 
                currentstructure = header
            elif state is None:
                print("error parsing RDfile ")

        if line.startswith('$DTYPE'):
            dtype = line[12:]

        if dtype and line.startswith('$DATUM'):
            datum = line[7:]

            if dtype.startswith('RX'):
                data[dtype] = datum
            elif dtype.startswith('TRANSFORM'):
                data[dtype[14:]] = datum
            dtype = None

    data['rxnsmiles'] = tosmiles(products, reactants)
    return data

def tosmiles(products, reactants):
    try:
       smiles = ''
       for r in reactants:
           mol = Chem.MolFromMolBlock(r)
           smiles += Chem.MolToSmiles(mol) + '.'
       smiles = smiles[:-1]
       smiles += '>>'
       for p in products:
           mol = Chem.MolFromMolBlock(r)
           smiles += Chem.MolToSmiles(mol) + '.'
       smiles = smiles[:-1]
   
       rxn = Chem.rdChemReactions.ReactionFromSmarts(smiles)
       sm = Chem.rdChemReactions.ReactionToSmiles(rxn)
       return sm
    except:
       return None

def readrdfiles(fname):
    """ read all of the individual SDFiles from the concatenated SDFile """
    print('readrdfiless', fname)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    count = 0
    conn=psql.connect(user=dbname)
    with gzip.open(fname, 'rt') as file:
        # remove header lines
        file.readline()
        file.readline()
        # loop over concatenated files
        while True:
            rdrecord = readnextRDfile(file)
            if not rdrecord:
                break;
            if 'RX_ID' in rdrecord.keys():
                writedb(conn, rdrecord)
                count += 1
                if count % 10000 == 0:
                   conn.commit()

    conn.commit()
    conn.close()
    print("wrote ", count, " records")




def writedb(conn, data):
     """ write a SDFile record the database """

     sql = 'insert into reaxys.rdfile (%s) values %s'
     with conn.cursor() as cur:
         columns = data.keys()
         values = [data[column] for column in columns]
         if debug:
             print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))))
         else:
             cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))



def readrdfile():
  """ read the SDFiles. This requires special functions because this is
    not an XML file
  """
  for filepath in glob.iglob('rdf/*.rdf.gz'):
    readrdfiles(filepath)

readrdfile()
