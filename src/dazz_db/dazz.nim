# vim: sw=2 ts=2 sts=2 tw=80 et:
from db import nil

proc Open*(path: string): ref db.DAZZ_DB =
  new(result)
  var modpath = path # var, not let, so Open_DB can alter it
  discard db.Open_DB(modpath, addr result[])

proc Close*(db1: ref db.DAZZ_DB) =
  db.Close_DB(addr db1[])
