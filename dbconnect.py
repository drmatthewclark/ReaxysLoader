import psycopg2 as psql
import os
from subprocess import check_output
from credentials import *  # import dbase, user, password, host variables

def getConnection():
    conn = psql.connect(dbname=dbase, user=user, password=password, host=host)
    return conn

def psql_cmd(command):
    env = dict(os.environ, PGPASSWORD=password)
    out = check_output(['psql', '-h', host, '-U', user, '-d', dbase, '-c', command], env=env)
    print('output:', out.decode('utf8'))

