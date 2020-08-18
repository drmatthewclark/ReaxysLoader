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

CHUNKSIZE = 10000

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

# end of loop over lines, process them
# data is a dictionary with elements like
    reacts, prods, data['rxnsmiles'] = tosmiles(products, reactants)
    
    sql = 'insert into reaxys.molecule (%s) values %s;'
    for i, (regno, smiles) in enumerate(reacts):
        name = ''
        for j in range(0 , len(data['RX_RXRN'])):
            rxn = data['RX_RXRN'][j]
            if rxn == regno and 'RX_RCT' in data.keys():
                name = data['RX_RCT'][j]
                break

        dbdata = {'molecule_id' : regno, 'name' : name, 'molstructure' : smiles}
        writerecord(conn, sql, dbdata)

    for i, (regno, smiles) in enumerate(prods):
        name = ''
        for j in range(0 , len(data['RX_PXRN'])):
            rxn = data['RX_PXRN'][j]
            if rxn == regno and 'RX_PRO' in data.keys():
                name = data['RX_PRO'][j]
                break

        dbdata = {'molecule_id' : regno, 'name' : name, 'molstructure' : smiles}
        writerecord(conn, sql, dbdata)

    if 'RX_RCT' in data.keys():
        del data['RX_RCT']
    if 'RX_PRO' in data.keys():
        del data['RX_PRO']
    #print("\twrote %7i molecule records" % (moleculecount))
    return data



def tosmiles(products, reactants):
    """ create smiles strings for products, reactants, and the reaction """
    """ pruducts, reactatnts are lists of tuples of regno,molfile """
    prods = list()
    reacts = list()
    smiles = ''
    for regno,r in reactants:
        try:
            mol = Chem.MolFromMolBlock(r, sanitize=False)
            reactant_smiles = (regno, Chem.CanonSmiles(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)))
        except:
            with open(efilename, 'a') as file:
                file.write("molecule\n%s\n error %s line %i\n" % ( r, sys.exc_info()[1], sys.exc_info()[2].tb_lineno)  )
            reactant_smiles = ''

        if reactant_smiles != '':
            reacts.append(reactant_smiles)
            smiles += reactant_smiles[1] + '.'

    smiles = smiles[:-1]
    smiles += '>>'

    for regno, p in products:
        try:
            mol = Chem.MolFromMolBlock(p, sanitize=False)
            product_smiles = (regno, Chem.CanonSmiles(Chem.MolToSmiles(mol, isomericSmiles = True, canonical = True)))
        except:
            with open(efilename, 'a') as file:
                file.write("molecule\n%s\n error %s line %i\n" % ( p, sys.exc_info()[1], sys.exc_info()[2].tb_lineno)  )
            product_smiles = ''

        if product_smiles != '':
            prods.append(product_smiles)
            smiles += product_smiles[1]  + '.'
            smiles = smiles[:-1]
    try:   
        rxn = Chem.rdChemReactions.ReactionFromSmarts(smiles)
        sm = Chem.rdChemReactions.ReactionToSmiles(rxn)
    except:
        sm = None
        with open(efilename, 'a') as file:
            file.write("molecule\n%s\n error %s line %i\n" % ( smiles, sys.exc_info()[1], sys.exc_info()[2].tb_lineno)  )

    return reacts, prods, sm


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

        with conn.cursor() as cur:
            cur.execute( '\n'.join(insertcache))

        insertcache.clear()
        conn.commit()

    print("\twrote %7i reaction records" %(count))



def writerecord(conn, sql, data):
     """ write a SDFile record the database """
     global hashset
     global insertcache

     count = 0
     with conn.cursor() as cur:
         columns = data.keys()
         values = [data[column] for column in columns]
         cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8')
         h = hash(cmd)

         if not h in hashset:
            hashset.add(h) 
            insertcache.add(cmd)
            count += 1

            if len(insertcache) > CHUNKSIZE:
                cur.execute( '\n'.join(insertcache))
                insertcache.clear()
                conn.commit()

     return count


            
def readrdfiles():
  """ read the SDFiles. This requires special functions because this is
    not an XML file
  """
  global insertcache
  oldlen = 0

  conn=psql.connect(user=dbname)

  for i, filepath in enumerate(glob.iglob('rdf/*.rdf.gz')):
        readrdfile(filepath, conn)
        newlen = len(hashset)
        new = newlen - oldlen
        print("\tadded %6i total records: %6i" %( new, newlen))
        oldlen = len(hashset)
  
  conn.commit()
  conn.close()


readrdfiles()
