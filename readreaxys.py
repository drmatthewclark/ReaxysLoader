# rmc file release number
rmcversion = 'rx200171'
dbname='mclark'
debug = False

import xml.etree.ElementTree as ET
import psycopg2 as psql
from   psycopg2 import sql
from psycopg2.extensions import AsIs
import glob
import gzip
import hashlib

# for rdkit Mol to Smiles
#from rdkit import Chem
#from rdkit import RDLogger

def getid(element):
    """ returns a hash code for the XML element and children to create a repeatable id for the element """
    if element is None:
        return -1

    id = element.attrib.get('ID') 
    if id and id.isnumeric() and len(list(element)) == 1:
       return int(id)
    else:
       text = ''.join(element.itertext()).encode('utf-8')
       hashcode = hash(hashlib.md5(text).hexdigest())
       if debug:
           print('key:',text, 'hash', hashcode)
       return hashcode

def readconditions(fname, dbname):
    """ 
    read an xml file into the designated database

    fname - xml file to read
    key - field for the object key
    dbname - name for the database to store data to
    sql - template sql statement to execute for storing the data
    """
    print('readconditions', fname)
    tree = ET.parse(gzip.open(fname));
    root = tree.getroot()
    conn = psql.connect(user=dbname)
    cur = conn.cursor()
    sql =  'insert into reaxys.conditions (%s) values %s on conflict (condition_id) do nothing';
    for record in root.findall('REACTIONS/REACTION'):
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
                minmax = {}
                if subelem:
                    for subsubelem in subelem:
                        tag = subsubelem.tag
                        value = subsubelem.text
                        if value == '-INF':
                          value = '-99999999'
                        elif value == 'INF':
                          value = '99999999'
                        minmax[tag] = value
                
                    data[condition] = '[' + minmax['min'] + ',' + minmax['max'] + ']'
 
            columns = data.keys()
            values = [data[column] for column in columns]
            if debug: 
                print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))))
            else:
                cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))
    conn.commit()
    conn.close()

def readstages(fname, dbname):
    """ 
    read an xml file into the designated database

    fname - xml file to read
    key - field for the object key
    dbname - name for the database to store data to
    sql - template sql statement to execute for storing the data
    """
    print('readstages', fname)
    tree = ET.parse(gzip.open(fname));
    root = tree.getroot()
    conn = psql.connect(user=dbname)
    cur = conn.cursor()
    sql =  'insert into reaxys.stages (%s) values %s on conflict (stage_id) do nothing';
    for record in root.findall('REACTIONS/REACTION'):
        for elem in record.findall('.//STAGES'): 
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
            if debug:
                print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))))
            else:
                cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))
    conn.commit()
    conn.close()

def readvariations(fname, dbname):
    """ 
    read an xml file into the designated database

    fname - xml file to read
    key - field for the object key
    dbname - name for the database to store data to
    sql - template sql statement to execute for storing the data
    """
    print('readvariations',fname)
    tree = ET.parse(gzip.open(fname));
    root = tree.getroot()
    conn = psql.connect(user=dbname)
    cur = conn.cursor()
    sql =  'insert into reaxys.variation (%s) values %s on conflict (variation_id) do nothing';
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
 
            for item  in ['STAGES', 'COMMENTS', 'IDENTIFIERS', 'LINKS', 'GROUPS', 'ANIMALS', 'CIT_ID', 'CONDITIONS', 'REACTANTS', 
                  'PRODUCTS', 'REAGENTS','CATALYSTS','SOLVENTS','METABOLITES','KEYWORDS']:
                subelem  = elem.findall(item)
                if subelem:
                    ilist = list()
                    for sselem in subelem:
                        tag = getid(sselem)
                        ilist.append(tag)
                    data[item] = ilist 

            columns = data.keys()
            values = [data[column] for column in columns]
            if debug:
                print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values)))) 
            else:
                cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))
    conn.commit()
    conn.close()

def readreactions(fname, dbname):
    """ 
    read an xml file into the designated database

    fname - xml file to read
    key - field for the object key
    dbname - name for the database to store data to
    sql - template sql statement to execute for storing the data
    """
    print('readreactions ',fname)
    tree = ET.parse(gzip.open(fname));
    root = tree.getroot()
    conn = psql.connect(user=dbname)
    cur = conn.cursor()
    sql =  'insert into reaxys.reaction (%s) values %s on conflict (reaction_id) do nothing';
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
                       text = hash(hashlib.md5(text.encode('utf-8')).hexdigest())
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
        if debug:
            print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values)))) 
        else:
            cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))
    conn.commit()
    conn.close()


def readsubstances(fname, dbname):
    """ 
    read an xml file into the designated database

    fname - xml file to read
    key - field for the object key
    dbname - name for the database to store data to
    sql - template sql statement to execute for storing the data
    """
    print('readsubstances', fname)
    tree = ET.parse(gzip.open(fname));
    root = tree.getroot()
    conn = psql.connect(user=dbname)
    cur = conn.cursor()
    sql =  'insert into reaxys.substance (%s) values %s on CONFLICT (substance_id) do nothing';
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
                    if debug:
                        print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))))
                    else:
                        cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))
    conn.commit()
    conn.close() 


def readcitations(fname, dbname):
    """ 
    read an xml file into the designated database

    fname - xml file to read
    key - field for the object key
    dbname - name for the database to store data to
    sql - template sql statement to execute for storing the data
    """
    print('readcitation', fname)
    tree = ET.parse(gzip.open(fname));
    root = tree.getroot()
    conn = psql.connect(user=dbname)
    cur = conn.cursor()
    sql =  'insert into reaxys.citation (%s) values %s on conflict(citation_id) do nothing';

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
        if debug:
            print(cur.mogrify(sql, (AsIs(','.join(columns)), tuple(values))))
        else:
            cur.execute(sql, (AsIs(','.join(columns)), tuple(values)))
    conn.commit()
    conn.close()


#for filepath in glob.iglob('udm-cit/*citations*.xml.gz'):
#  readcitations(filepath, 'mclark')

go = False
for filepath in glob.iglob('udm-rea/*reactions*.xml.gz'):
        readreactions(filepath, 'mclark')
        readconditions(filepath, 'mclark')
        readstages(filepath, 'mclark')
        readvariations(filepath, 'mclark')
        readsubstances(filepath, 'mclark')
