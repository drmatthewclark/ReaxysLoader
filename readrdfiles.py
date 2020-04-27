# rmc file release number
dbname='mclark'
debug = False

import xml.etree.ElementTree as ET
import psycopg2 as psql
from psycopg2.extensions import AsIs
import glob
import gzip

# for rdkit Mol to Smiles
from rdkit import Chem
from rdkit import RDLogger
import rdkit.Chem.rdChemReactions

global conn


def readnextRDfile(file):
    """  read the next reaction from the concatenated file """
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
    """ process the separated RDfile into components, and smiles """
    global counter
    # header to make molfile parser happy
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
            dtype = line[line.find('RX'):]

        if dtype and line.startswith('$DATUM'):
            datum = line[7:]

            if dtype.endswith('RCT') or dtype.endswith('PRO'):
                if not dtype in data.keys():
                    data[dtype] = list()
                data[dtype].append(datum) 
            elif dtype.endswith('XRN'):
                if not dtype in data.keys():
                    data[dtype] = list()
                data[dtype].append(int(datum)) 
            elif dtype.startswith('RX'):
                data[dtype] = datum
            elif dtype.startswith('TRANSFORM'):
                data[dtype[14:]] = datum
            dtype = None

    reacts, prods, data['rxnsmiles'] = tosmiles(products, reactants)
    
    sql = 'insert into reaxys.molecule (molecule_id, name, molstructure) values (%s, %s, %s) on conflict (molecule_id) do nothing'

    if data and 'RX_RXRN' in data.keys():
        for i in range(0 , len(data['RX_RXRN'])):
            rid = data['RX_RXRN'][i] 
            smiles = None
            name = None
            if reacts is not None and len(reacts) > i:
                 smiles = reacts[i]
            if 'RX_RCT' in data.keys() and len(data['RX_RCT']) > i:
                name = data['RX_RCT'][i]

            with conn.cursor() as cur:
                if debug:
                    print(cur.mogrify(sql, (rid, name, smiles)))
                cur.execute(sql, (rid, name, smiles))

    if data and 'RX_PXRN' in data.keys():
        for i in range(0 , len(data['RX_PXRN'])):
            rid = data['RX_PXRN'][i]
            smiles = None
            name = None
            if prods is not None and len(prods) > i:
                 smiles = prods[i]
            if 'RX_PRO' in data.keys() and len(data['RX_PRO']) > i:
                 name = data['RX_PRO'][i]

            with conn.cursor() as cur:
                if debug:
                    print(cur.mogrify(sql, (rid, name, smiles)))
                cur.execute(sql, (rid, name, smiles))

    if 'RX_RCT' in data.keys():
        del data['RX_RCT']
    if 'RX_PRO' in data.keys():
        del data['RX_PRO']

    return data



def tosmiles(products, reactants):
    """ create smiles strings for products, reactants, and the reaction """
    prods = list()
    reacts = list()
    try:
       smiles = ''
       for r in reactants:
           mol = Chem.MolFromMolBlock(r)
           reactant_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical = True)
           reacts.append(reactant_smiles)
           smiles += reactant_smiles + '.'

       smiles = smiles[:-1]
       smiles += '>>'

       for p in products:
           mol = Chem.MolFromMolBlock(p)
           product_smiles = Chem.MolToSmiles(mol, isomericSmiles = True, canonical = True)
           prods.append(product_smiles)
           smiles += product_smiles + '.'
       smiles = smiles[:-1]
   
       rxn = Chem.rdChemReactions.ReactionFromSmarts(smiles)
       sm = Chem.rdChemReactions.ReactionToSmiles(rxn)
       return reacts, prods, sm
    except:
       return None, None, None



def readrdfiles(fname):
    """ read all of the individual SDFiles from the concatenated SDFile """
    print('readrdfiless', fname)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    count = 0
    global conn
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
                if count % 10000000 == 0:
                   conn.commit()

    conn.commit()
    conn.close()
    print("wrote ", count, " records")




def writedb(conn, data):
     """ write a SDFile record the database """

     sql = 'insert into reaxys.rdfile (%s) values %s on conflict (rx_id) do nothing'
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

