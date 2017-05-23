import DB

{.compile: "../../repos/DAZZ_DB/DB.c".}
{.compile: "../../repos/DAZZ_DB/QV.c".}
{.passC: "-I../../repos/DAZZ_DB".}


when isMainModule:
  echo "hi"
