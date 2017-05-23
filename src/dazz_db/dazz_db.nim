import db

{.compile: "../../repos/DAZZ_DB/DB.c".}
{.compile: "../../repos/DAZZ_DB/QV.c".}
{.passC: "-I../../repos/DAZZ_DB".}


when isMainModule:
  from os import execShellCmd
  var mydb: db.HITS_DB = db.HITS_DB()
  block:
    var badpath = "/badpath"
    doAssert -1 == db.Open_DB(badpath, addr mydb)
  var goodpath = "tmpdb.db"
  block:
    let ret = os.execShellCmd("fasta2DB " & goodpath & " -ifoo < /dev/null")
    doAssert 0 == ret
  block:
    let ret = db.Open_DB(goodpath, addr mydb)
    db.Close_DB(addr mydb)
    doAssert 0 == ret
  echo "DAZZ_DB test passed!"
