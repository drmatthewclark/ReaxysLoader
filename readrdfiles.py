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
# for rdkit Mol to Smiles
from rdkit import Chem
from rdkit import RDLogger
import rdkit.Chem.rdChemReactions
import hashlib

global hashset
global insertcache
index = 0

CHUNKSIZE = 20000
TEST_MODE = True

efilename='structures.errors'
hashset = set()
insertcache = set()
lastline = None

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
    
    if line.startswith('$RFMT $RIREG'):
        rdfile = line
        line = file.readline()
    else:
        rdfile = lastline

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

    lastline = line
    tags = processRXN(rdfile, conn)
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
    dtype = None
    data = {}
    regno = 0

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
            prefix = line[12:rx-1] # RCT(2) to track number and name
            #$DTYPE ROOT:RCT(2):RX_RXRN

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
    
    sql = 'insert into reaxys.molecule (%s) values %s;'
    for i, (regno, smiles) in enumerate(reacts):
        name = ''
        for j in range(0 , len(data['RX_RXRN'])):
            prefix, rxn = data['RX_RXRN'][j]
            if rxn == regno and 'RX_RCT' in data.keys():
                names = data['RX_RCT']
                for j in range(0, len(names)):
                    pfix, tname = names[j]
                    if ( pfix == prefix):
                         name = tname
                break
       
        (regno,unpacked) = reactants[i]
             
        dbdata = {'molecule_id' : regno, 'name' : name, 'smiles' : smiles, 'sdfile' : unpacked }
 
        writerecord(conn, sql, dbdata)

    for i, (regno, smiles) in enumerate(prods):
        name = ''
        for j in range(0 , len(data['RX_PXRN'])):
            prefix, rxn = data['RX_PXRN'][j]
            if rxn == regno and 'RX_PRO' in data.keys():
                names = data['RX_PRO']
                for j in range(0, len(names)):
                        pfix, tname  = names[j]
                        if (pfix == prefix):
                            name = tname
                break

        (regno,unpacked) = products[i]

        dbdata = {'molecule_id' : regno, 'name' : name, 'smiles' : smiles, 'sdfile' : unpacked }
        writerecord(conn, sql, dbdata)

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

    #print("\twrote %7i molecule records" % (moleculecount))

    return data

# centralize this function
def createSmiles(molblock):

    smiles = None
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

        if reactant_smiles != (regno, None):
            reacts.append(reactant_smiles)
            smiles += reactant_smiles[1] + '.'

    smiles = smiles[:-1]
    smiles += '>>'

    for regno, p in products:
        product_smiles = (regno, createSmiles(p)) 

        if product_smiles != (regno, None):
            prods.append(product_smiles)
            smiles += product_smiles[1]  + '.'

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
        sql = 'insert into reaxys.rdfile (%s) values %s;'
        while True:
            rdrecord = readnextRDfile(file, conn)
            if not rdrecord:
                break;
            if 'RX_ID' in rdrecord.keys():
                count += writerecord(conn, sql, rdrecord)
        
        if not TEST_MODE:
            with conn.cursor() as cur:
                cur.execute( '\n'.join(insertcache))
        
        insertcache.clear()
        conn.commit()

    print("\twrote %7i reaction records" %(count))


"""
  wrapper for hash to remove sdfile from hashing 
  so that it is not considered for uniquifying records

"""
def hashrecord(record):

    if type(record) is dict and 'sdfile' in record.keys():
        temp = record.copy() # python is object oriented
        del temp['sdfile']
        return hash(temp)

    return hash(record)



def writerecord(conn, sql, data):
     """ write a SDFile record the database """
     global hashset
     global insertcache

     count = 0
     with conn.cursor() as cur:
         columns = data.keys()
         values = [data[column] for column in columns]
         cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8')
         if TEST_MODE:
            print("--start record")
            print(cmd)
            print("--end record")

         h = hashrecord(cmd)

         if not h in hashset:
            hashset.add(h) 
            insertcache.add(cmd)
            count += 1
            if len(insertcache) > CHUNKSIZE:
                if not TEST_MODE:
                    cur.execute( '\n'.join(insertcache))

                insertcache.clear()
                conn.commit()

     return count


            
def readrdfiles():
  """ read the SDFiles. This requires special functions because this is
    not an XML file
  """
  global insertcache
  global index
  conn=psql.connect(user=dbname)
  key = "_"

  for i, filepath in enumerate(glob.iglob('rdf/*.rdf.gz')):
        index = os.path.basename(filepath)
        index = int(index[index.find(key) + len(key):-7 ])
        print('file index', index)
        oldlen = len(hashset)
        readrdfile(filepath, conn)
        newlen = len(hashset)
        new = newlen - oldlen
        print("\tadded %6i total records: %6i" %( new, newlen))
  
  conn.commit()
  conn.close()


readrdfiles()
