# ReaxysLoader
Load Reaxys fileset into a postgresql database

* now corrected hash function for ID creation. This can be tested by adding primary key after table loading. It now does not complain about duplicate primary keys.

create postgresql database<br>
create schemas with psql < loader_schema (creates reaxys schema)<br> 
load database:<br>
python readreaxys.py  - loads all but structure files<br>
python readrdfiles.py - loads structure files as reaction smiles<br>
<br>
add indices:  psql < loader_index  (adding after load makes loading faster)

UDM data structures are created as objects, with an MD5 hash key created on the contents<br>

