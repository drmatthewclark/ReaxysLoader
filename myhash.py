import hashlib

def myhash(text):
    """ replaces the built in python hash function with one that is deterministic """
 
    # use md5 to create a hash
    digest  = hashlib.md5(text).hexdigest()
    dint = int(digest,16)  # create 128 bit integer

    # fold hi/low bytes into a 64 bit integer.
    hibytes = ((dint & 0xFFFFFFFFFFFFFFFF0000000000000000) >> 64)
    lobytes =   dint & 0x0000000000000000FFFFFFFFFFFFFFFF
    # xor high and low halves
    result  = (hibytes ^ lobytes) 
    # if greater than max 64bit signed int, make negative. python doesn't wrap on overflow
    # because it never overflows
    if (result > 0x7FFFFFFFFFFFFFFF):
        result = -(result >> 1) 

    return result
