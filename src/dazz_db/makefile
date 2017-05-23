default:
	nim c -r DB.nim
QV.nim:
	c2nim --header --cdecl --out:$@ ../../repos/DAZZ_DB/QV.h
DB.nim:
	c2nim --header --cdecl --out:$@ ../../repos/DAZZ_DB/DB.h
#DB.h:
#	cpp ../../repos/DAZZ_DB/DB.h > DB.h
#	# Then modify.

# You still have to remove include-guard constant and un-indent.
