# vim: sw=2 ts=2 sts=2 tw=80 et:

from dazz_db/db import nil

#discard db.Malloc(1024, "Allocated?")
var mydb: db.HITS_DB
let name = "foo"
discard db.Open_DB(name, addr mydb) # can over-write name
echo "Opened ", name
