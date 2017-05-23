for x in ../../repos/DAZZ_DB/*.c; do echo "{.compile: \"$(basename $x)\".}"; done

