#
# dedup some substance id's not previously deduplictated and append varying names
#
from dbconnect import getConnection
import psycopg2 as psql

conn = getConnection()

dups   =  "select substance_id  from (select count(substance_id) as count, substance_id from reaxys.substance group by substance_id)a where a.count >1;"
getdup =  'select ctid,* from reaxys.substance where substance_id = %s order by name asc;'
update =  'update reaxys.substance set name = %s where ctid = %s;'
delete =  'delete from reaxys.substance where ctid = %s;'

cur = conn.cursor()
cur.execute(dups)
for line in cur.fetchall():
    with conn.cursor() as c2:
        c2.execute(getdup, line)
        main=c2.fetchone()
        firstname=str(main[3])
        firstctid=main[0]

        for l in c2.fetchall():
            newname = str(l[3])
            ctid=l[0]
            if newname != 'None':
                firstname=firstname+ '|' + str(l[3])
            
            print(c2.mogrify(update, (firstname, firstctid)))
            print(c2.mogrify(delete, (ctid,)))
        
            c2.execute(update, (firstname, firstctid))
            c2.execute(delete, (ctid,))

conn.commit()
conn.close()
