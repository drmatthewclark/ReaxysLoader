import hashlib

def myhash(text):
    """ replaces the built in python hash function with one that is deterministic """
    digest  = hashlib.md5(text).hexdigest()
    dint = int(digest,16)  # create 128 bit integer
    # use the first bit as a sign bit
    # fold hi/low bytes into a 64 bit integer.
    hibytes = ((dint & 0xFFFFFFFFFFFFFFFF0000000000000000) >> 64)
    lobytes =   dint & 0x0000000000000000FFFFFFFFFFFFFFFF
    # xor high and low halves
    result  = (hibytes ^ lobytes) 
    # if greater than max 64bit signed int, make negative. python doesn't wrap
    if (result > 0x7FFFFFFFFFFFFFFF):
        result = -1 *  (result >> 1) 

    return result
