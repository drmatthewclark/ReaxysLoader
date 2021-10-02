#!/usr/bin/python
#  run the processes

import os
from dbconnect import getConnection

mydir = os.path.dirname(os.path.realpath(__file__))

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
                print(command)
                cur.execute(command)
                conn.commit()

    conn.close()



import readreaxys
import readrdfiles
import substance_dedup
sqlfile(mydir + '/reaxys_index')

print('reaxys load complete')
