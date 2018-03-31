# vim: sw=2 ts=2 sts=2 tw=80 et:

from dazz_db/db import nil
from dazz_db/dazz import nil

proc main() =
  #discard db.Malloc(1024, "Allocated?")
  var db1: db.DAZZ_DB
  let name = "foo"
  discard db.Open_DB(name, addr db1) # can over-write name
  echo "Opened ", name
  db.Close_DB(addr db1)

  let db2 = dazz.Open(name)
  defer: dazz.Close(db2)
  echo "Re-Opened ", name

main()
