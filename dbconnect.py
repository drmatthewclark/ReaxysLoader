import psycopg2 as psql


dbname='mclark'

def getConnection():
    conn=psql.connect(user=dbname)
    return conn


