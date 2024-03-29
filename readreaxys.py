# database to load data to
# if true will print out the sql statements and other data for debugging
debug = False

import xml.etree.ElementTree as ET
import psycopg2 as psql
from psycopg2.extensions import AsIs
import glob
import tempfile
import os
import time
import gzip
from myhash import myhash
from pathlib import Path
from dbconnect import getConnection

path = Path('.')
version = os.path.basename(path.parent.absolute())
print('loading version', version)
mydir = os.path.dirname(os.path.realpath(__file__))


CHUNKSIZE = 50000

lines = set()

def getid(element):
    """ returns a hash code for the XML element and children to create a repeatable id for the element """

    if element is None:
        print("getid error: null element")
        return -1

    id = element.attrib.get('ID') 

    if id and id.isnumeric() and len(list(element)) == 1:
       return int(id)
    else:
       text = ET.tostring(element)
       hashcode = myhash(text)

       if debug:
           print('key:',text, 'hash', hashcode)
       return hashcode


def readconditions(tree, conn):
    """ 
    read an xml file into the designated database
    """
    start = time.time() 
    insertcache = set()
    root = tree.getroot()

    sql =  'insert into reaxys_temp.conditions (%s) values %s;'
    cur = conn.cursor()
    for record in root.findall('REACTIONS'):
        for elem in record.findall('.//CONDITIONS'):
            data = {}
            data['condition_id'] = getid(elem) 

            for condition in ['ATMOSPHERE','PREPARATION','REFLUX']:
                 subelem = elem.findall(condition)
                 if subelem:
                    for sselem in subelem:
                        tag = sselem.tag
                        value = sselem.text
                        data[tag] = value
    
            #this set has ranges
            for condition in ['PH','PRESSURE','REACTION_MOLARITY','TEMPERATURE','TIME','TOTAL_VOLUME']:
                subelem = elem.find(condition)
                if subelem:
                    mint = subelem.find('min')
                    maxt = subelem.find('max')
                    minval = mint.text
                    maxval = maxt.text
                    if minval == '-INF':
                        minval = '-1e100'
                    if maxval == 'INF':
                        maxval = '1e100'
            
                    array =  '[' + minval + ',' + maxval + ']'
                    data[condition] = array
 
            columns = data.keys()
            values = [data[column] for column in columns]
            cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8') + "\n"

            h = hash(cmd)
            if not h in lines:
                lines.add(h)
                insertcache.add(cmd)
                if len(insertcache) > CHUNKSIZE:
                   cur.execute( '\n'.join(insertcache))
                   insertcache.clear() 

    if len(insertcache) > 0:
        cur.execute( '\n'.join(insertcache))

    cur.close()
    conn.commit()
    print("\treadconditions load took %5.2f %6i records" % ((time.time() - start), len(lines)))



def readstages(tree, conn):
    """ 
    read an xml file into the designated database
    """
    start = time.time() 
    insertcache = set()
    root = tree.getroot()

    sql =  'insert into reaxys_temp.stages (%s) values %s;'
    cur = conn.cursor()
    for record in root.findall('REACTIONS/REACTION'):
        for elem in record.findall('VARIATIONS/STAGES'): 
            data = {}
            data['stage_id'] = getid(elem) 
    
            #this set has ranges
            for item  in ['CONDITIONS', 'REACTANTS', 'PRODUCTS', 'REAGENTS','CATALYSTS','SOLVENTS','METABOLITES']:
                subelem  = elem.findall(item)
                if subelem:
                    ilist = list()
                    for sselem in subelem:
                        tag = getid(sselem)
                        ilist.append(tag)
                    data[item] = ilist
 
            columns = data.keys()
            values = [data[column] for column in columns]
            cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8') + "\n"  

            h = hash(cmd)
            if not h in lines:
                lines.add(h)
                insertcache.add(cmd)
                if len(insertcache) > CHUNKSIZE:
                   cur.execute( '\n'.join(insertcache))
                   insertcache.clear() 

    if len(insertcache) > 0:
        cur.execute( '\n'.join(insertcache))

    cur.close()
    conn.commit()
    print("\treadstages load took %5.2f %6i records" % ((time.time() - start), len(lines)))



def readvariations(tree, conn):
    """ 
    read an xml file into the designated database
    """
    start = time.time() 
    insertcache = set()
    root = tree.getroot()

    sql =  'insert into reaxys_temp.variation (%s) values %s;'
    cur = conn.cursor()
    for record in root.findall('REACTIONS/REACTION'):
        for elem in record.findall('VARIATIONS'): 
            data = {}
            data['variation_id'] = getid(elem) 

            #this set has lists
            for item  in ['CREATION_DATE', 'EXPERIMENT_ID','EXPERIMENT_TYPE','MODIFICATION_DATE','PROJECT_NAME','QUALIFICATION',
                'SOURCE','DESTINATION','CONCLUSION_PHRASE','CREATED_AT_SITE','DUPLICATE_EXP_REF','PREPARATIVE','ELN_CITATION',
                'REACTION_SCALE','NEXTMOVE_REACTION_TYPE','RXNO_REACTION_TYPE','ANALYTICAL_DATA_EXISTS']:
                subelem  = elem.find(item)
                if subelem is not None:
                    tag = subelem.tag
                    value = subelem.text
                    data[tag] = value
 
            for item  in ['STAGES', 'COMMENTS', 'IDENTIFIERS', 'LINKS', 'GROUPS', 'ANIMALS', 'CONDITIONS', 'REACTANTS', 
                  'PRODUCTS', 'REAGENTS','CATALYSTS','SOLVENTS','METABOLITES','KEYWORDS']:
                subelem  = elem.findall(item)
                if subelem:
                    ilist = list()
                    for sselem in subelem:
                        tag = getid(sselem)
                        ilist.append(tag)
                    data[item] = ilist 

            for item  in [ 'CIT_ID' ]:
                subelem  = elem.findall(item)
                if subelem:
                    ilist = list()
                    for sselem in subelem:
                        tag = int(sselem.text)
                        ilist.append(tag)
                    data[item] = ilist 

            columns = data.keys()
            values = [data[column] for column in columns]
            cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8') + "\n"

            h = hash(cmd)
            if not h in lines:
                lines.add(h)
                insertcache.add(cmd)
                if len(insertcache) > CHUNKSIZE:
                   cur.execute( '\n'.join(insertcache))
                   insertcache.clear() 

    if len(insertcache) > 0:
        cur.execute( '\n'.join(insertcache))

    cur.close()
    conn.commit()
    print("\treadvariations load took %5.2f %6i records" % ((time.time() - start), len(lines)))
    return



def readreactions(tree, conn):
    """ 
    read an xml file into the designated database
    """
    start = time.time() 
    insertcache = set()
    root = tree.getroot()

    sql =  'insert into reaxys_temp.reaction (%s) values %s;'
    cur = conn.cursor()

    for elem in root.findall('REACTIONS/REACTION'):
        data = {}
        data['reaction_id'] = getid(elem) 
        data['reaxys_reaction_id'] = elem.attrib.get('ID');

        for item in ['RXNSTRUCTURE','RANK','MW_LARGEST_PRODUCT','SORT_CREATION_DATE','SORT_REACTION_SCALE','SORT_TOTAL_VOLUME']:
            subelem  = elem.find(item)
            if subelem is not None:
                tag = subelem.tag
                value = subelem.text
                data[tag] = value
 
        for item in ['REACTANT_ID', 'METABOLITE_ID', 'PRODUCT_ID']:
            subelem  = elem.findall(item)
            if subelem:
                ilist = list()
                for sselem in subelem:
                    # some of the data are MD5 hash strings, to be annoying
                    # we will just make this a hash to be compatible as integer key 
                    # for this database
                    text  = sselem.text
                    if text.startswith("MD5"):
                       text = myhash(text.encode('utf-8'))
                    text = int(text)
                    ilist.append(text)
                data[item] = ilist 

        for item in ['VARIATIONS']:
            subelem  = elem.findall(item)
            if subelem:
                ilist = list()
                for sselem in subelem:
                    tag = getid(sselem)
                    ilist.append(tag)
                data[item] = ilist 

        columns = data.keys()
        values = [data[column] for column in columns]
        cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8') + "\n"

        h = hash(cmd)
        if not h in lines:
           lines.add(h)
           insertcache.add(cmd)
           if len(insertcache) > CHUNKSIZE:
              cur.execute( '\n'.join(insertcache))
              insertcache.clear() 

    if len(insertcache) > 0:
        cur.execute( '\n'.join(insertcache))

    cur.close()
    conn.commit()
    print("\treadreactions load took %5.2f %6i records" % ((time.time() - start), len(lines)))
    return




def readsubstances(tree, conn):
    """ 
    read an xml file into the designated database
    """
    start = time.time() 
    insertcache = set()
    root = tree.getroot()

    sql =  'insert into reaxys_temp.substance (%s) values %s;'
    cur = conn.cursor()

    for record in root.findall('REACTIONS/REACTION'):
        for elem in record.findall('VARIATIONS'): 

            #this set has lists
            for item  in ['REACTANTS', 'PRODUCTS', 'REAGENTS','CATALYSTS','SOLVENTS','METABOLITES']:
                data = {}
                subelem  = elem.find(item)
                if subelem:
                    data['substance_id'] = getid(subelem)
                    temp = subelem.attrib.get('ID')
                    if temp.isnumeric():
                       data['reaxys_id'] = temp
                       
                    for subsubelem in subelem:
                        tag = subsubelem.tag
                        value = subsubelem.text
                        data[tag] = value 

                    columns = data.keys()
                    values = [data[column] for column in columns]
                    cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode('utf-8') + "\n" 

                    h = hash(cmd)
                    if not h in lines:
                        lines.add(h)
                        insertcache.add(cmd)
                        if len(insertcache) > CHUNKSIZE:
                            cur.execute( '\n'.join(insertcache))
                            insertcache.clear() 

    if len(insertcache) > 0:
        cur.execute( '\n'.join(insertcache))

    cur.close()
    conn.commit()
    print("\treadsubstances load took %5.2f %6i records" % ((time.time() - start), len(lines)) )
    return


def readcitations(tree, conn):
    """ 
    read an xml file into the designated database
    """
    insertcache = set()
    start = time.time() 

    root = tree.getroot()
 
    sql =  'insert into reaxys_temp.citation (%s) values %s;'
    cur = conn.cursor()

    for elem in root.findall('CITATIONS/CITATION'):
        data = {}
        id = elem.attrib.get('ID')
        data['citation_id'] = id 
        
        for subelem in elem:
            value = subelem.text
            tag = subelem.tag
            data[tag] = value

        columns = data.keys()
        values = [data[column] for column in columns]

        cmd = cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))).decode() + '\n'
        h = hash(cmd)
        if not h in lines:
           lines.add(h)
           insertcache.add(cmd)
           if len(insertcache) > CHUNKSIZE:
              cur.execute( '\n'.join(insertcache))
              insertcache.clear() 

    if len(insertcache) > 0:
        cur.execute( '\n'.join(insertcache))

    cur.close()
    conn.commit()
    print("\treadcitations load took %5.2f %6i records" % ((time.time() - start), len(lines)))
    return



def sqlfile(fname):
    """ read and execute a sql file"""
    c = getConnection()
    with open(fname, 'r') as f:
        sql = f.read()

    commands = sql.split(';')

    print('executing sql', fname)
    with c.cursor() as cur:
        for command in commands:
            command = command.strip()
            if command != '':
                cur.execute(sql)
                c.commit()
    c.close()


def load():
    
    conn = getConnection()
    sqlfile(mydir + '/reaxys_schema')

    with conn.cursor() as cur:
        cur.execute('insert into reaxys_temp.version (version) values (%s);', (version,))
        conn.commit()
 
    file_list = glob.iglob('udm-cit/*citations*.xml.gz')
    file_list.sort()

    for i, filepath in enumerate(file_list):
        print("file: ", filepath)
        tree = ET.parse(gzip.open(filepath));
        readcitations(tree, conn)

    lines.clear()

    file_list = glob.iglob('udm-rea/*reactions*.xml.gz')
    file_list.sort()
     
    for i, filepath in enumerate(file_list):
       print("file: ", filepath)
       tree = ET.parse(gzip.open(filepath));
       readreactions(tree, conn)
       readconditions(tree, conn)
       readstages(tree, conn)
       readvariations(tree, conn)
       readsubstances(tree, conn)

load()
