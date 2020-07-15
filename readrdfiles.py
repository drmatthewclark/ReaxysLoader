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
from dedup import writedb

global SFILE 
global hashset
SFILE = None
filename='structures.table2'
efilename='structures.errors'
hashset = set()

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
                # add molfile to list
                reactants.append(currentstructure)
                currentstructure = header
            elif state == 'product':
                currentstructure += tail
                # add molfile to list
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

# end of loop over lines, process them
    reacts, prods, data['rxnsmiles'] = tosmiles(products, reactants)
    
    sql = 'insert into reaxys.molecule (%s) values %s;'

    if data and reacts and 'RX_RXRN' in data.keys():
        numreagents= len(data['RX_RXRN'])
        assert numreagents >= len(reacts)
        for i in range(0 , len(reacts)):
            rid = data['RX_RXRN'][i] 
            smiles = None
            name = None
            if reacts is not None and len(reacts) > i:
                 smiles = reacts[i]
            if 'RX_RCT' in data.keys() and len(data['RX_RCT']) > i:
                name = data['RX_RCT'][i]
            if smiles is not None:
                dbdata = {'molecule_id' : rid, 'name' : name, 'molstructure' : smiles}
                writerecord(conn, sql, dbdata)

    if data and prods and 'RX_PXRN' in data.keys():
        assert len(data['RX_PXRN']) >= len(prods)
        for i in range(0 , len(data['RX_PXRN'])):
            rid = data['RX_PXRN'][i]
            smiles = None
            name = None
            if prods is not None and len(prods) > i:
                 smiles = prods[i]
            if 'RX_PRO' in data.keys() and len(data['RX_PRO']) > i:
                 name = data['RX_PRO'][i]
            if smiles is not None:
                dbdata = {'molecule_id' : rid, 'name' : name, 'molstructure' : smiles}
                writerecord(conn, sql, dbdata)

    if 'RX_RCT' in data.keys():
        del data['RX_RCT']
    if 'RX_PRO' in data.keys():
        del data['RX_PRO']
    #print("\twrote %7i molecule records" % (moleculecount))
    return data



def tosmiles(products, reactants):
    """ create smiles strings for products, reactants, and the reaction """
    prods = list()
    reacts = list()
    smiles = ''
    for r in reactants:
        try:
            mol = Chem.MolFromMolBlock(r, sanitize=False)
            reactant_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        except:
            with open(efilename, 'a') as file:
                file.write("molecule\n%s\n error %s line %i\n" % ( r, sys.exc_info()[1], sys.exc_info()[2].tb_lineno)  )
            reactant_smiles = ''

        reacts.append(reactant_smiles)
        smiles += reactant_smiles + '.'

    smiles = smiles[:-1]
    smiles += '>>'

    for p in products:
        try:
            mol = Chem.MolFromMolBlock(p, sanitize=False)
            product_smiles = Chem.MolToSmiles(mol, isomericSmiles = True, canonical = True)
        except:
            with open(efilename, 'a') as file:
                file.write("molecule\n%s\n error %s line %i\n" % ( p, sys.exc_info()[1], sys.exc_info()[2].tb_lineno)  )
            product_smiles = ''

        prods.append(product_smiles)
        smiles += product_smiles + '.'
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

    print("\twrote %7i reaction records" %(count))



def writerecord(conn, sql, data):
     """ write a SDFile record the database """
     global SFILE
     global hashset
     count = 0
     with conn.cursor() as cur:
         columns = data.keys()
         values = [data[column] for column in columns]
         cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values)))
         h = hash(cmd)

         if not h in hashset:
            SFILE.write(cmd.decode('utf-8')+'\n')
            hashset.add(h) 
            count += 1

     return count


            
def readrdfiles():
  """ read the SDFiles. This requires special functions because this is
    not an XML file
  """
  global SFILE
  SFILE =  open(filename,"w")   
  oldlen = 0

  conn=psql.connect(user=dbname)

  for i, filepath in enumerate(glob.iglob('rdf/*.rdf.gz')):
        readrdfile(filepath, conn)
        newlen = len(hashset)
        new = newlen - oldlen
        print("\tadded %6i total records: %6i" %( new, newlen))
        oldlen = len(hashset)

  SFILE.close()

  writedb(filename, conn)
  conn.close()

readrdfiles()
