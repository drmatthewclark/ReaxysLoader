import psycopg2 as psql


dbname='xmclark'

def getConnection():
    conn=psql.connect(user=dbname)
    return conn


