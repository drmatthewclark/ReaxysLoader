# ReaxysLoader
Load Reaxys fileset into a postgresql database

create postgresql database<br>
create schemas with psql < schema; psql < rdfileschema<br>
load database:<br>
python readreaxys.py  - all but structure files<br>
python readrdfiles.py - loads structure files as reaction smiles<br>
<br>
UDM data structures are created as objects, with an MD5 hash key created on the contents<br>

