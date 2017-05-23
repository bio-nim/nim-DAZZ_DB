# vim: sw=2 ts=2 sts=2 tw=80 et:
import dazz_db/dazz_db
import dazz_db/db

from os import execShellCmd

echo "DAZZ_DB test running ..."

var mydb: db.HITS_DB = db.HITS_DB()
block:
  var badpath = "/badpath"
  doAssert -1 == db.Open_DB(badpath, addr mydb)
let goodpath = "tmpdb.db"
block:
  doAssert 0 == os.execShellCmd("fasta2DB " & goodpath & " -ifoo < /dev/null")
block:
  doAssert 0 == db.Open_DB(goodpath, addr mydb)
  db.Close_DB(addr mydb)
block:
  doAssert 0 == os.execShellCmd("DBrm " & goodpath)

echo "DAZZ_DB test passed!"
