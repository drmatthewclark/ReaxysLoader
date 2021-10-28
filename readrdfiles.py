# rmc file release number
dbname='mclark'
debug = False
import xml.etree.ElementTree as ET
import psycopg2 as psql
from psycopg2.extensions import AsIs
import glob
import gzip
import os
import sys
import time
from dbconnect import getConnection

# for rdkit Mol to Smiles
from rdkit import Chem
from rdkit import RDLogger
import rdkit.Chem.rdChemReactions
import hashlib

global hashset
global insertcache
index = 0

CHUNKSIZE = 50000
TEST_MODE = False
STORE_RDFILE = False

hashset = set()
insertcache = set()
lastline = None
counts = {}
counts['molecule'] = 0
counts['reaction'] = 0
counts['other'] = 0

rireg = len('$RFMT $RIREG')


def readnextRDfile(file, conn):
    """  read the next reaction from the concatenated file """
    line = '' # file.readline()
    rdfile = ''
    tags = {}
    blankcount = 0
    global lastline

    line = file.readline().rstrip()

    # these lines only at beginning of file
    if (line.startswith('$RDFILE') or line.startswith('$DATM')):
        line == readline()

    if line == '':
        return None
    
    if line[:rireg] =='$RFMT $RIREG':
        rdfile = line
        line = file.readline()
    else:
        rdfile = lastline

    while not line[:rireg] == '$RFMT $RIREG':
        line = file.readline()
        if not line[:rireg] == '$RFMT $RIREG':
            rdfile = rdfile + line
        # python file reading is brain dead and returns an empty string
        # at the EOF instad of  None or any other EOF marker so we have to
        # guess that a bunch of blank lines means that we have reached the
        # end of file 
        if line == '':
            blankcount += 1
        if blankcount > 2:
            break

    lastline = line
    tags = processRXN(rdfile, conn)
    # add the whole RDFILE to the database
    # it is quite large and not terribly useful for most purposes
    if STORE_RDFILE:
        tags['rx_rdfile'] = rdfile[:-1] #trim last char
    tags['rx_file_id'] = index
    return tags # dictionary 


global counter 
counter = 0



def processRXN(rdfile, conn):
    """ process the separated RDfile into components, and smiles """
    global counter
    # header to make molfile parser happy
    header = '\n' + 'GSMACCS-II07189510252D 1   0.00366     0.00000     0' + '\n\n'  + '  0  0  0     0  0            999 V3000' + '\n'
    tail = 'M  END' 
    counter += 1
    reactants = list()
    products = list()
    currentstructure = header
    state = ''
    name = '' 
    dtype = None
    data = {}
    regno = 0
    msql = 'insert into reaxys_temp.molecule (%s) values %s;'

    for line in rdfile.splitlines():
        endctab = None
        # check this first for performance
        if line.startswith("M  "):
            if line.startswith('M  V30 BEGIN RE'): # reactant section
                state = 'reactant'  # the structures are reactants
                regno = 0
                continue
            elif line.startswith('M  V30 END RE'):
                state = None   # reactant section has ended
                continue
            elif line.startswith('M  V30 BEGIN PR'):
                state = 'product'  # molecules are now products
                regno = 0
                continue
            elif line.startswith('M  V30 END PR'):
                state = None  # product section has ended
                continue
            elif line.startswith('M  V30 END CT'):
                endctab = True 
            elif line.startswith('M  V30 CO'):
                equals = line.find('=')
                if equals != -1:
                    #registy no. of current molecule
                    regno = int(line[equals+1:].strip())
                else:
                    regno = 0


# example:
#   M  V30 COUNTS 14 14 0 0 0 REGNO=747939
        if state:
            currentstructure += line + '\n'
       
        if endctab:
            endctab = None
            if state == 'reactant':
                currentstructure += tail
                # add molfile to list
                reactants.append((regno,currentstructure))
                currentstructure = header
            elif state == 'product':
                currentstructure += tail
                # add molfile to list
                products.append((regno,currentstructure))
                currentstructure = header
            elif state is None:
                print("error parsing RDfile ")

        if line.startswith('$DTYPE'):
            rx = line.find('RX')
            dtype = line[rx:]
            prefix = line[:rx] 

        if dtype and line.startswith('$DATUM'):
            datum = line[7:]

            if dtype.endswith('RCT') or dtype.endswith('PRO'):
                if not dtype in data.keys():
                    data[dtype] = list()
                data[dtype].append((prefix,datum))  # name
            elif dtype.endswith('XRN'):
                if not dtype in data.keys():
                    data[dtype] = list()
                data[dtype].append((prefix, int(datum))) # reg number
            elif dtype.startswith('RX'):
                data[dtype] = datum
            elif dtype.startswith('TRANSFORM'):
                data[dtype[14:]] = datum
            dtype = None

# end of loop over lines, process them
# data is a dictionary with elements like
    reacts, prods, data['rxnsmiles'] = tosmiles(products, reactants)

    for i, (regno, smiles) in enumerate(reacts):
        for (rxprefix, rxn) in data['RX_RXRN']:
            if rxn == regno and 'RX_RCT' in data.keys():
                name = ''
                for (rctprefix, tname) in data['RX_RCT']:
                    if rctprefix == rxprefix:
                        name = tname
                        break
 
        (regno,unpacked) = reactants[i]
        dbdata = {'molecule_id' : regno, 'name' : name, 'smiles' : smiles, 'sdfile' : unpacked, 'rx_file_id': index}
 
        writerecord(conn, msql, dbdata)

 
    for i, (regno, smiles) in enumerate(prods):
        for (rxprefix, rxn) in data['RX_PXRN']:
            if rxn == regno and 'RX_PRO' in data.keys():
                name = ''
                for (rctprefix, tname) in data['RX_PRO']:
                    if rctprefix == rxprefix:
                        name = tname
                        break

        (regno,unpacked) = products[i]

        dbdata = {'molecule_id' : regno, 'name' : name, 'smiles' : smiles, 'sdfile' : unpacked, 'rx_file_id': index}
        writerecord(conn, msql, dbdata)

    # remove the names, those are in the molecule table
    if 'RX_RCT' in data.keys():
        del data['RX_RCT']
    if 'RX_PRO' in data.keys():
        del data['RX_PRO']

    # switch tuples to registry numbers
    for i in range(0, len(data['RX_PXRN'])):
        prefix, rxn = data['RX_PXRN'][i]
        data['RX_PXRN'][i] = rxn
        
    for i in range(0, len(data['RX_RXRN'])):
        prefix, rxn = data['RX_RXRN'][i]
        data['RX_RXRN'][i] = rxn

    return data

# centralize this function
def createSmiles(molblock):

    smiles = '' 

    if molblock is None:
        return smiles
    try:
        mol = Chem.MolFromMolBlock(molblock, strictParsing=False, sanitize=True)
        if mol is not None:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    except:
        pass

    return smiles


def tosmiles(products, reactants):
    """ create smiles strings for products, reactants, and the reaction """
    """ pruducts, reactatnts are lists of tuples of regno,molfile """
    prods = list()
    reacts = list()
    smiles = ''

    for regno,r in reactants:
        reactant_smiles = (regno, createSmiles(r)) 
        reacts.append(reactant_smiles)
        if reactant_smiles[1] != '':
            smiles += reactant_smiles[1] + '.'

    smiles = smiles[:-1]
    smiles += '>>'

    for regno, p in products:
        product_smiles = (regno, createSmiles(p)) 
        prods.append(product_smiles)
        if product_smiles[1] != '':
            smiles += product_smiles[1] + '.'

    smiles = smiles[:-1]

    return reacts, prods, smiles


def readrdfile(fname, conn):

    """ read all of the individual SDFiles from the concatenated SDFile """
    print('readrdfiles: ', fname)
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    count = 0
    with gzip.open(fname, 'rt') as file:
        # remove header lines
        file.readline()
        file.readline()
        # loop over concatenated files
        sql = 'insert into reaxys_temp.rdfile (%s) values %s;'
        while True:
            rdrecord = readnextRDfile(file, conn)
            if not rdrecord:
                break;
            if 'RX_ID' in rdrecord.keys():
                count += writerecord(conn, sql, rdrecord)

    flush(conn) 

    print("\tprocessed %7i reaction records" %(count))



def flush(conn):
        if not TEST_MODE and len(insertcache) > 0:
            with conn.cursor() as cur:
                cur.execute( '\n'.join(insertcache))
        insertcache.clear()
        conn.commit()


"""
  wrapper for hash to remove sdfile from hashing 
  so that it is not considered for uniquifying records

"""
def hashrecord(record):

    if 'sdfile' in record:
        temp = record.split("\n")[0]
        return hash(temp)

    return hash(record)



def writerecord(conn, sql, data):
     """ write a SDFile record the database """
     global hashset
     global insertcache
     rectype = None     
     count = 0

     with conn.cursor() as cur:
         columns = data.keys()
         values = [data[column] for column in columns]
         cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8')
         if TEST_MODE:
            print("--start record")
            print(cmd)
            print("--end record")
         
         if 'molecule_id' in columns:
            hashdata = 'm' + str(data['molecule_id'])
            rectype = 'molecule'
         elif 'RX_ID' in columns:
            hashdata = 'r' + str(data['RX_ID'])
            rectype = 'reaction'
         else:
            hashdata = 'o'+ cmd
            rectype = 'other'

         h = hashrecord(hashdata)

         if not h in hashset:
            hashset.add(h) 
            insertcache.add(cmd)
            count += 1
            counts[rectype] += 1
            if len(insertcache) > CHUNKSIZE:
                flush(conn)

     return count


def sqlfile(fname):
    """ read and execute a sql file"""
    conn = getConnection()
    with open(fname, 'r') as f:
        sql = f.read()

    commands = sql.split(';')
    
    print('executing sql', fname)
    with conn.cursor() as cur:
        for command in commands:
            command = command.strip()
            if command != '':
                cur.execute(sql)
                conn.commit()

    conn.close()

            
def readrdfiles():
  """ read the SDFiles. This requires special functions because this is
    not an XML file
  """
  global insertcache
  global index

  conn = getConnection()

  key = "_"
  numfiles = len(glob.glob('rdf/*.rdf.gz'))

  for i, filepath in enumerate(glob.iglob('rdf/*.rdf.gz')):
        start = time.time()
        index = os.path.basename(filepath)
        index = int(index[index.find(key) + len(key):-7 ])
        print('file index', index)
        oldlen = len(hashset)
        readrdfile(filepath, conn)
        newlen = len(hashset)
        new = newlen - oldlen
        elapsed = time.time() - start 
        remaining = '%5.2f' % ((numfiles-i-1)*elapsed )
        elapsed = "%5.2f" %(elapsed)
        print('\ttook:',elapsed, 'remaining:', remaining, 'inserts:',counts)
  
  conn.commit()
  conn.close()

readrdfiles()
