default:
	nim c --debuginfo -r dazz.nim
QV.nim:
	c2nim --header --cdecl --out:$@ ../../repos/DAZZ_DB/QV.h
DB.nim:
	c2nim --header --cdecl --out:$@ --skipinclude ../../repos/DAZZ_DB/DB.h
#DB.h:
#	cpp ../../repos/DAZZ_DB/DB.h > DB.h
#	# Then modify.

# You still have to remove include-guard constant and un-indent.

clean:
	rm -rf nimcache/
foo: CFLAGS=-I../../repos/DAZZ_DB
foo: CPPFLAGS=-g
foo: foo.c ../../repos/DAZZ_DB/DB.c ../../repos/DAZZ_DB/QV.c
