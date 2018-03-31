VPATH=../repos/DAZZ_DB
PATH:=../repos/DAZZ_DB:${PATH}
export PATH

default: foo.db
	${MAKE} test-open
test-open:
	nim c -r --debuginfo --debugger:native -o:open.exe open
foo.db: DBrm
	-DBrm -v foo.db
	fasta2DB -v foo.db -ifoo < foo.fasta
	ls -larth .foo*
DBrm:
	cd ../repos/DAZZ_DB; make -j
	which DBrm
	which fasta2DB
clean:
	rm -f foo.db .foo.* open.exe
