import os
import xml.etree.ElementTree as ET

def dedup(fname):
    outfilename = fname + "_tmp" # tempfile.mkstemp()[1]
    lines_seen = set() # holds lines already seen
    outfile = open(outfilename, "w")
    removed = 0
    written = 0
    for line in open(fname, "r"):
        h = hash(line)    
        if h not in lines_seen: # not a duplicate
         outfile.write(line)
         lines_seen.add(h)
         written += 1
        else:
         removed += 1

    outfile.close()
    os.remove(fname)
    os.rename(outfilename, fname)
    print("%s removed %d lines wrote %d" % (fname, removed, written))
    return outfilename


#
# write data from the file to the database
#
def writedb(fname, conn):

    cur = conn.cursor()
    blocksize = 50000
    dedup(fname)
    statement = ''
    count = 0
    for line in open(fname, 'r'):
        count += 1
        statement += line
        if (count % blocksize == 0):
            cur.execute(statement)
            conn.commit()
            statement = ''

    cur.execute(statement)
    cur.close()
    conn.commit()


