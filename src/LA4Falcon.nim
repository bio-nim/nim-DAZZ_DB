##  vim: set et ts=2 sts=2 sw=2 :
## 
##   Utility for displaying the overlaps in a .las file in a variety of ways including
##     a minimal listing of intervals, a cartoon, and a full out alignment.
## 
##   Author:    Gene Myers
##   Creation:  July 2013
##   Last Mod:  Jan 2015
## 
##   Based on the original LAshow.c, this code is modified by Jason Chin to support generating
##     consensus sequences from daligner output

import
  DB, DBX, align, license_myers

const
  MAX_OVERLAPS* = 50000

#template MIN*(X, Y: untyped): untyped =
#  if ((X) < (Y)): (X) else: (Y)

var GROUP*: bool = false

##  Allows us to group overlaps between a pair of a/b reads as a unit, one per
##  direction (if applicable).  beg/end will point to the same overlap when
##  only one overlap found.

type
  OverlapGroup* = object
    beg*: Overlap
    `end`*: Overlap
    score*: cint
    blen*: cint


var ovlgrps*: ptr OverlapGroup

proc compare_ovlgrps*(grp1: pointer; grp2: pointer): cint =
  return (cast[ptr OverlapGroup](grp2)).score -
      (cast[ptr OverlapGroup](grp1)).score

proc belongs*(grp: ptr OverlapGroup; ovl: ptr Overlap): bool =
  var prev: ptr Overlap = addr(grp.`end`)
  return prev.flags == ovl.flags and (ovl.path.abpos > prev.path.aepos) and
      (ovl.path.bbpos > prev.path.bepos) and
      (ovl.path.abpos - prev.path.aepos) < 251

##  Add a new overlap to a new or existing overlap group. Always adds when group
##  flag is false, effectively greating groups of 1.
##  Returns 1 if added as a new overlap group, otherwise 0.
##  caller keeps track of count

proc add_overlap*(aln: ptr Alignment; ovl: ptr Overlap; count: cint): bool =
  var added: cint = false
  ##  we assume breads are in order
  if not GROUP or count < 0 or ovlgrps[count].beg.bread != ovl.bread:
    ##  Haven't seen this bread yet (or we're not grouping), move to new overlap group
    var next: ptr OverlapGroup = addr(ovlgrps[count + 1])
    next.beg = ovl[]
    next.`end` = ovl[]
    next.blen = aln.blen
    var p: ptr Path = addr(ovl.path)
    var olen: cint = p.bepos - p.bbpos
    var hlen: cint = (MIN(p.abpos, p.bbpos)) +
        (MIN(aln.alen - p.aepos, aln.blen - p.bepos))
    next.score = olen - hlen
    added = true
  else:
    var curr: ptr OverlapGroup = addr(ovlgrps[count])
    ##  Seen, should we combine it with the previous overlap group or move
    ##  on to the next?
    if belongs(curr, ovl):
      curr.`end` = ovl[]
      ##  rescore
      var beg: ptr Overlap = addr(curr.beg)
      var `end`: ptr Overlap = addr(curr.`end`)
      var olen: cint = `end`.path.bepos - beg.path.bbpos
      var hlen: cint = (MIN(beg.path.abpos, beg.path.bbpos)) +
          (MIN(aln.alen - `end`.path.aepos, aln.blen - `end`.path.bepos))
      curr.score = olen - hlen
    else:
      var next: ptr OverlapGroup = addr(ovlgrps[count + 1])
      next.beg = ovl[]
      next.`end` = ovl[]
      next.blen = aln.blen
      var p: ptr Path = addr(ovl.path)
      var olen: cint = p.bepos - p.bbpos
      var hlen: cint = (MIN(p.abpos, p.bbpos)) +
          (MIN(aln.alen - p.aepos, aln.blen - p.bepos))
      next.score = olen - hlen
      added = true
  return added

proc print_hits*(hit_count: cint; dbx2: ptr HITS_DBX; bbuffer: cstring;
                buffer: ptr char; bsize: int64; MAX_HIT_COUNT: cint) =
  var tmp_idx: cint
  qsort(ovlgrps, (hit_count + 1), sizeof((OverlapGroup)), compare_ovlgrps)
  tmp_idx = 0
  while tmp_idx < (hit_count + 1) and tmp_idx < MAX_HIT_COUNT:
    var grp: ptr OverlapGroup = addr(ovlgrps[tmp_idx])
    ## Load_ReadX assuming db2 == db1 is true
    Load_ReadX(dbx2, grp.`end`.bread, bbuffer, 0)
    if COMP(grp.`end`.flags): Complement_Seq(bbuffer, grp.blen)
    Upper_Read(bbuffer)
    var rlen: int64 = (int64)(grp.`end`.path.bepos) - (int64)(grp.beg.path.bbpos)
    if rlen < bsize:
      strncpy(buffer, bbuffer + grp.beg.path.bbpos, rlen)
      buffer[rlen - 1] = '\0'
      printf("%08d %s\x0A", grp.`end`.bread, buffer)
    else:
      fprintf(stderr,
              "[WARNING]Skipping super-long read %08d, len=%lld, buf=%lld\x0A",
              grp.`end`.bread, rlen, bsize)
    inc(tmp_idx)
  printf("+ +\x0A")

var Usage*: ptr cstring = ["[-mfsocargUFM] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ", "    <src1:db|dam> [ <src2:db|dam> ] <align:las> [ <reads:FILE> | <reads:range> ... ]"]

const
  LAST_READ_SYMBOL* = '$'

proc ORDER*(l: pointer; r: pointer): cint =
  var x: cint = (cast[ptr int32](l))[]
  var y: cint = (cast[ptr int32](r))[]
  return x - y

proc main*(argc: cint; argv: ptr cstring): cint =
  var
    _dbx1: HITS_DBX
    dbx1: ptr HITS_DBX = addr(_dbx1)
  var
    _dbx2: HITS_DBX
    dbx2: ptr HITS_DBX = addr(_dbx2)
  var db1: ptr HITS_DB = addr(dbx1.db)
  var db2: ptr HITS_DB = addr(dbx2.db)
  var
    _ovl: Overlap
    ovl: ptr Overlap = addr(_ovl)
  var
    _aln: Alignment
    aln: ptr Alignment = addr(_aln)
  var input: ptr FILE
  var novl: int64
  var
    tspace: cint
    tbytes: cint
    small: cint
  var
    reps: cint
    pts: ptr cint
  var input_pts: cint
  var
    ALIGN: cint
    CARTOON: cint
    REFERENCE: cint
    FLIP: cint
  var
    INDENT: cint
    WIDTH: cint
    BORDER: cint
    UPPERCASE: cint
  var ISTWO: cint
  var MAP: cint
  var
    FALCON: cint
    OVERLAP: cint
    M4OVL: cint
  ##  XXX: MAX_HIT_COUNT should be renamed
  var
    SEED_MIN: cint
    MAX_HIT_COUNT: cint
    SKIP: cint
  var PRELOAD: cint
  ##   Process options
  var
    i: cint
    j: cint
    k: cint
  var flags: array[128, cint]
  var eptr: cstring
  ARG_INIT("LA4Falcon")
  INDENT = 4
  WIDTH = 100
  BORDER = 10
  FALCON = 0
  M4OVL = 0
  SEED_MIN = 8000
  SKIP = 0
  ALIGN = 0
  REFERENCE = 0
  CARTOON = 0
  FLIP = 0
  MAX_HIT_COUNT = 400
  j = 1
  i = 1
  while i < argc:
    if argv[i][0] == '-':
      case argv[i][1]
      of 'i':
        ARG_NON_NEGATIVE(INDENT, "Indent")
      of 'w':
        ARG_POSITIVE(WIDTH, "Alignment width")
      of 'b':
        ARG_NON_NEGATIVE(BORDER, "Alignment border")
      of 'H':
        ARG_POSITIVE(SEED_MIN, "seed threshold (in bp)")
      of 'n':
        ARG_POSITIVE(MAX_HIT_COUNT, "max numer of supporting read ouput (used for FALCON consensus. default 400, max: 2000)")
        if MAX_HIT_COUNT > 2000: MAX_HIT_COUNT = 2000
      else:
        ARG_FLAGS("smfocargUFMP")
    else:
      argv[inc(j)] = argv[i]
    inc(i)
  argc = j
  UPPERCASE = flags['U']
  ALIGN = flags['a']
  REFERENCE = flags['r']
  CARTOON = flags['c']
  FLIP = flags['F']
  MAP = flags['M']
  OVERLAP = flags['o']
  M4OVL = flags['m']
  FALCON = flags['f']
  SKIP = flags['s']
  GROUP = flags['g']
  PRELOAD = flags['P']
  ##  Preload DB reads, if possible.
  if argc <= 2:
    fprintf(stderr, "Usage: %s %s\x0A", Prog_Name, Usage[0])
    fprintf(stderr, "       %*s %s\x0A", cast[cint](strlen(Prog_Name)), "", Usage[1])
    exit(1)
  ##   Open trimmed DB or DB pair
  var status: cint
  var
    pwd: cstring
    root: cstring
  var input: ptr FILE
  ISTWO = 0
  status = Open_DBX(argv[1], dbx1, PRELOAD)
  if status < 0: exit(1)
  if db1.part > 0:
    fprintf(stderr, "%s: Cannot be called on a block: %s\x0A", Prog_Name, argv[1])
    exit(1)
  if argc > 3:
    pwd = PathTo(argv[3])
    root = Root(argv[3], ".las")
    if (input = fopen(Catenate(pwd, "/", root, ".las"), "r")) != nil:
      ISTWO = 1
      fclose(input)
      status = Open_DBX(argv[2], dbx2, PRELOAD)
      if status < 0: exit(1)
      if db2.part > 0:
        fprintf(stderr, "%s: Cannot be called on a block: %s\x0A", Prog_Name,
                argv[2])
        exit(1)
      Trim_DB(db2)
    else:
      dbx2 = dbx1
      db2 = db1
    free(root)
    free(pwd)
  else:
    dbx2 = dbx1
    db2 = db1
  Trim_DB(db1)
  ##   Process read index arguments into a sorted list of read ranges
  input_pts = 0
  if argc == ISTWO + 4:
    if argv[ISTWO + 3][0] != LAST_READ_SYMBOL or argv[ISTWO + 3][1] != '\0':
      var
        eptr: cstring
        fptr: cstring
      var
        b: cint
        e: cint
      b = strtol(argv[ISTWO + 3], addr(eptr), 10)
      if eptr > argv[ISTWO + 3] and b > 0:
        if eptr[] == '-':
          if eptr[1] != LAST_READ_SYMBOL or eptr[2] != '\0':
            e = strtol(eptr + 1, addr(fptr), 10)
            input_pts = (fptr <= eptr + 1 or fptr[] != '\0' or e <= 0)
        else:
          input_pts = (eptr[] != '\0')
      else:
        input_pts = 1
  if input_pts:
    var
      v: cint
      x: cint
    var input: ptr FILE
    input = Fopen(argv[ISTWO + 3], "r")
    if input == nil: exit(1)
    reps = 0
    while (v = fscanf(input, " %d", addr(x))) != EOF:
      if v == 0:
        fprintf(stderr, "%s: %d\'th item of input file %s is not an integer\x0A",
                Prog_Name, reps + 1, argv[2])
        exit(1)
      else:
        inc(reps, 1)
    reps = reps * 2
    pts = cast[ptr cint](Malloc(sizeof((int) * reps), "Allocating read parameters"))
    if pts == nil: exit(1)
    rewind(input)
    v = 0
    while v < reps:
      fscanf(input, " %d", addr(x))
      pts[v] = pts[v + 1] = x
      inc(v, 2)
    fclose(input)
  else:
    pts = cast[ptr cint](Malloc(sizeof((int) * 2 * argc), "Allocating read parameters"))
    if pts == nil: exit(1)
    reps = 0
    if argc > 3 + ISTWO:
      var
        c: cint
        b: cint
        e: cint
      var
        eptr: cstring
        fptr: cstring
      c = 3 + ISTWO
      while c < argc:
        if argv[c][0] == LAST_READ_SYMBOL:
          b = db1.nreads
          eptr = argv[c] + 1
        else:
          b = strtol(argv[c], addr(eptr), 10)
        if eptr > argv[c]:
          if b <= 0:
            fprintf(stderr, "%s: %d is not a valid index\x0A", Prog_Name, b)
            exit(1)
          if eptr[] == '\0':
            pts[inc(reps)] = b
            pts[inc(reps)] = b
            continue
          elif eptr[] == '-':
            if eptr[1] == LAST_READ_SYMBOL:
              e = INT32_MAX
              fptr = eptr + 2
            else:
              e = strtol(eptr + 1, addr(fptr), 10)
            if fptr > eptr + 1 and fptr[] == 0 and e > 0:
              pts[inc(reps)] = b
              pts[inc(reps)] = e
              if b > e:
                fprintf(stderr, "%s: Empty range \'%s\'\x0A", Prog_Name, argv[c])
                exit(1)
              continue
        fprintf(stderr, "%s: argument \'%s\' is not an integer range\x0A",
                Prog_Name, argv[c])
        exit(1)
        inc(c)
      qsort(pts, reps div 2, sizeof((int64)), ORDER)
      b = 0
      c = 0
      while c < reps:
        if b > 0 and pts[b - 1] >= pts[c] - 1:
          if pts[c + 1] > pts[b - 1]: pts[b - 1] = pts[c + 1]
        else:
          pts[inc(b)] = pts[c]
          pts[inc(b)] = pts[c + 1]
        inc(c, 2)
      pts[inc(b)] = INT32_MAX
      reps = b
    else:
      pts[inc(reps)] = 1
      pts[inc(reps)] = INT32_MAX
  ##   Initiate file reading and read (novl, tspace) header
  var
    over: cstring
    pwd: cstring
    root: cstring
  pwd = PathTo(argv[2 + ISTWO])
  root = Root(argv[2 + ISTWO], ".las")
  over = Catenate(pwd, "/", root, ".las")
  input = Fopen(over, "r")
  if input == nil: exit(1)
  if fread(addr(novl), sizeof((int64)), 1, input) != 1: SYSTEM_ERROR
  if fread(addr(tspace), sizeof((int)), 1, input) != 1: SYSTEM_ERROR
  if tspace == 0:
    printf("\x0ACRITICAL ERROR: tspace=0 in \'%s\'", root)
    exit(1)
  if tspace <= TRACE_XOVR:
    small = 1
    tbytes = sizeof((uint8))
  else:
    small = 0
    tbytes = sizeof((uint16))
  if not (FALCON or M4OVL):
    printf("\x0A%s: ", root)
    Print_Number(novl, 0, stdout)
    printf(" records\x0A")
  free(pwd)
  free(root)
  ##   Read the file and display selected records
  var j: cint
  var trace: ptr uint16
  var work: ptr Work_Data
  var tmax: cint
  var
    `in`: cint
    npt: cint
    idx: cint
    ar: cint
  var tps: int64
  var p_aread: int64 = - 1
  var buffer: array[131072, char]
  var skip_rest: cint = 0
  var
    abuffer: cstring
    bbuffer: cstring
  var
    ar_wide: cint
    br_wide: cint
  var
    ai_wide: cint
    bi_wide: cint
  var
    mn_wide: cint
    mx_wide: cint
  var tp_wide: cint
  var
    blast: cint
    match: cint
    seen: cint
    lhalf: cint
    rhalf: cint
  var hit_count: cint
  aln.path = addr((ovl.path))
  if ALIGN or REFERENCE or FALCON:
    work = New_Work_Data()
    abuffer = New_Read_Buffer(db1)
    bbuffer = New_Read_Buffer(db2)
    if FALCON:
      ovlgrps = calloc(sizeof((OverlapGroup)), MAX_OVERLAPS + 1)
      hit_count = - 1
  else:
    abuffer = nil
    bbuffer = nil
    work = nil
  tmax = 1000
  trace = cast[ptr uint16](Malloc(sizeof((uint16) * tmax), "Allocating trace vector"))
  if trace == nil: exit(1)
  `in` = 0
  npt = pts[0]
  idx = 1
  ar_wide = Number_Digits(cast[int64](db1.nreads))
  br_wide = Number_Digits(cast[int64](db2.nreads))
  ai_wide = Number_Digits(cast[int64](db1.maxlen))
  bi_wide = Number_Digits(cast[int64](db2.maxlen))
  if db1.maxlen < db2.maxlen:
    mn_wide = ai_wide
    mx_wide = bi_wide
    tp_wide = Number_Digits(cast[int64](db1.maxlen div tspace) + 2)
  else:
    mn_wide = bi_wide
    mx_wide = ai_wide
    tp_wide = Number_Digits(cast[int64](db2.maxlen div tspace) + 2)
  inc(ar_wide, (ar_wide - 1) div 3)
  inc(br_wide, (br_wide - 1) div 3)
  inc(ai_wide, (ai_wide - 1) div 3)
  inc(bi_wide, (bi_wide - 1) div 3)
  inc(mn_wide, (mn_wide - 1) div 3)
  inc(tp_wide, (tp_wide - 1) div 3)
  if FLIP:
    var x: cint
    x = ar_wide
    ar_wide = br_wide
    br_wide = x
    x = ai_wide
    ai_wide = bi_wide
    bi_wide = x
  blast = - 1
  match = 0
  seen = 0
  lhalf = rhalf = 0
  j = 0
  while j < novl:
    Read_Overlap(input, ovl)
    if ovl.path.tlen > tmax:
      tmax = (cast[cint](1.2 * ovl.path.tlen)) + 100
      trace = cast[ptr uint16](Realloc(trace, sizeof((uint16) * tmax),
                                    "Allocating trace vector"))
      if trace == nil: exit(1)
    ovl.path.trace = cast[pointer](trace)
    Read_Trace(input, ovl, tbytes)
    ##   Determine if it should be displayed
    ar = ovl.aread + 1
    if `in`:
      while ar > npt:
        npt = pts[inc(idx)]
        if ar < npt:
          `in` = 0
          break
        npt = pts[inc(idx)]
    else:
      while ar >= npt:
        npt = pts[inc(idx)]
        if ar <= npt:
          `in` = 1
          break
        npt = pts[inc(idx)]
    if not `in`:
      continue
    aln.alen = db1.reads[ovl.aread].rlen
    aln.blen = db2.reads[ovl.bread].rlen
    aln.flags = ovl.flags
    tps = ((ovl.path.aepos - 1) div tspace - ovl.path.abpos div tspace)
    if OVERLAP and not FALCON:
      if ovl.path.abpos != 0 and ovl.path.bbpos != 0: continue
      if ovl.path.aepos != aln.alen and ovl.path.bepos != aln.blen: continue
    if MAP:
      while ovl.bread != blast:
        if not match and seen and not (lhalf and rhalf):
          printf("Missing ")
          Print_Number(cast[int64](blast) + 1, br_wide + 1, stdout)
          printf(" %d ->%lld\x0A", db2.reads[blast].rlen, db2.reads[blast].coff)
        match = 0
        seen = 0
        lhalf = rhalf = 0
        inc(blast, 1)
      seen = 1
      if ovl.path.abpos == 0: rhalf = 1
      if ovl.path.aepos == aln.alen: lhalf = 1
      if ovl.path.bbpos != 0 or ovl.path.bepos != aln.blen: continue
      match = 1
    if not (FALCON or M4OVL):
      if ALIGN or CARTOON or REFERENCE: printf("\x0A")
      if FLIP:
        Flip_Alignment(aln, 0)
        Print_Number(cast[int64](ovl.bread) + 1, ar_wide + 1, stdout)
        printf("  ")
        Print_Number(cast[int64](ovl.aread) + 1, br_wide + 1, stdout)
      else:
        Print_Number(cast[int64](ovl.aread) + 1, ar_wide + 1, stdout)
        printf("  ")
        Print_Number(cast[int64](ovl.bread) + 1, br_wide + 1, stdout)
      if COMP(ovl.flags): printf(" c")
      else: printf(" n")
      printf("   [")
      Print_Number(cast[int64](ovl.path.abpos), ai_wide, stdout)
      printf("..")
      Print_Number(cast[int64](ovl.path.aepos), ai_wide, stdout)
      printf("] x [")
      Print_Number(cast[int64](ovl.path.bbpos), bi_wide, stdout)
      printf("..")
      Print_Number(cast[int64](ovl.path.bepos), bi_wide, stdout)
      printf("]")
    if M4OVL:
      var
        bbpos: int64
        bepos: int64
      var acc: cdouble
      if COMP(ovl.flags):
        bbpos = cast[int64](aln.blen) - cast[int64](ovl.path.bepos)
        bepos = cast[int64](aln.blen) - cast[int64](ovl.path.bbpos)
      else:
        bbpos = cast[int64](ovl.path.bbpos)
        bepos = cast[int64](ovl.path.bepos)
      acc = 100 -
          (200.0 * ovl.path.diffs) div
          (ovl.path.aepos - ovl.path.abpos + ovl.path.bepos - ovl.path.bbpos)
      printf("%09lld %09lld %lld %5.2f ", cast[int64](ovl.aread),
             cast[int64](ovl.bread), cast[int64](bbpos) - cast[int64](bepos), acc)
      printf("0 %lld %lld %lld ", cast[int64](ovl.path.abpos),
             cast[int64](ovl.path.aepos), cast[int64](aln.alen))
      printf("%d %lld %lld %lld ", COMP(ovl.flags), bbpos, bepos,
             cast[int64](aln.blen))
      if (cast[int64](aln.blen) < cast[int64](aln.alen)) and
          (cast[int64](ovl.path.bbpos) < 1) and
          (cast[int64](aln.blen) - cast[int64](ovl.path.bepos) < 1):
        printf("contains\x0A")
      elif (cast[int64](aln.alen) < cast[int64](aln.blen)) and
          (cast[int64](ovl.path.abpos) < 1) and
          (cast[int64](aln.alen) - cast[int64](ovl.path.aepos) < 1):
        printf("contained\x0A")
      else:
        printf("overlap\x0A")
    if FALCON:
      if p_aread == - 1:
        Load_ReadX(dbx1, ovl.aread, abuffer, 2)
        printf("%08d %s\x0A", ovl.aread, abuffer)
        p_aread = ovl.aread
        skip_rest = 0
      if p_aread != ovl.aread:
        print_hits(hit_count, dbx2, bbuffer, buffer, cast[int64](sizeof((buffer))),
                   MAX_HIT_COUNT)
        hit_count = - 1
        Load_ReadX(dbx1, ovl.aread, abuffer, 2)
        printf("%08d %s\x0A", ovl.aread, abuffer)
        p_aread = ovl.aread
        skip_rest = 0
      if skip_rest == 0:
        if add_overlap(aln, ovl, hit_count): inc(hit_count)
        if (hit_count + 1) > MAX_OVERLAPS: skip_rest = 1
        if false:
          tps = ((ovl.path.aepos - 1) div tspace - ovl.path.abpos div tspace)
          if small: Decompress_TraceTo16(ovl)
          Load_ReadX(dbx1, ovl.aread, abuffer, 0)
          Load_ReadX(dbx2, ovl.bread, bbuffer, 0)
          if COMP(aln.flags): Complement_Seq(bbuffer, aln.blen)
          Compute_Trace_PTS(aln, work, tspace)
          var tlen: cint = aln.path.tlen
          var trace: ptr cint = aln.path.trace
          var u: cint
          printf(" ")
          u = 0
          while u < tlen:
            printf("%d,", cast[int16](trace[u]))
            inc(u)
        if SKIP == 1:
          ## if SKIP = 0, then skip_rest is always 0
          if (cast[int64](aln.alen) < cast[int64](aln.blen)) and
              (cast[int64](ovl.path.abpos) < 1) and
              (cast[int64](aln.alen) - cast[int64](ovl.path.aepos) < 1):
            printf("* *\x0A")
            skip_rest = 1
    if ALIGN or CARTOON or REFERENCE:
      if ALIGN or REFERENCE:
        var
          aseq: cstring
          bseq: cstring
        var
          amin: cint
          amax: cint
        var
          bmin: cint
          bmax: cint
        if FLIP: Flip_Alignment(aln, 0)
        if small: Decompress_TraceTo16(ovl)
        amin = ovl.path.abpos - BORDER
        if amin < 0: amin = 0
        amax = ovl.path.aepos + BORDER
        if amax > aln.alen: amax = aln.alen
        if COMP(aln.flags):
          bmin = (aln.blen - ovl.path.bepos) - BORDER
          if bmin < 0: bmin = 0
          bmax = (aln.blen - ovl.path.bbpos) + BORDER
          if bmax > aln.blen: bmax = aln.blen
        else:
          bmin = ovl.path.bbpos - BORDER
          if bmin < 0: bmin = 0
          bmax = ovl.path.bepos + BORDER
          if bmax > aln.blen: bmax = aln.blen
        aseq = Load_Subread(db1, ovl.aread, amin, amax, abuffer, 0)
        bseq = Load_Subread(db2, ovl.bread, bmin, bmax, bbuffer, 0)
        aln.aseq = aseq - amin
        if COMP(aln.flags):
          Complement_Seq(bseq, bmax - bmin)
          aln.bseq = bseq - (aln.blen - bmax)
        else:
          aln.bseq = bseq - bmin
        Compute_Trace_PTS(aln, work, tspace, GREEDIEST)
        if FLIP:
          if COMP(aln.flags):
            Complement_Seq(aseq, amax - amin)
            Complement_Seq(bseq, bmax - bmin)
            aln.aseq = aseq - (aln.alen - amax)
            aln.bseq = bseq - bmin
          Flip_Alignment(aln, 1)
      if CARTOON:
        printf("  (")
        Print_Number(tps, tp_wide, stdout)
        printf(" trace pts)\x0A\x0A")
        Alignment_Cartoon(stdout, aln, INDENT, mx_wide)
      else:
        printf(" :   = ")
        Print_Number(cast[int64](ovl.path.diffs), mn_wide, stdout)
        printf(" diffs  (")
        Print_Number(tps, tp_wide, stdout)
        printf(" trace pts)\x0A")
      if REFERENCE:
        Print_Reference(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide)
      if ALIGN:
        Print_Alignment(stdout, aln, work, INDENT, WIDTH, BORDER, UPPERCASE, mx_wide)
    elif not (FALCON or M4OVL):
      printf(" :   < ")
      Print_Number(cast[int64](ovl.path.diffs), mn_wide, stdout)
      printf(" diffs  (")
      Print_Number(tps, tp_wide, stdout)
      printf(" trace pts)\x0A")
    inc(j)                    ##   Read it in
  if FALCON and hit_count != - 1:
    print_hits(hit_count, dbx2, bbuffer, buffer, cast[int64](sizeof((buffer))),
               MAX_HIT_COUNT)
    printf("- -\x0A")
    free(ovlgrps)
  free(trace)
  if ALIGN or FALCON:
    free(bbuffer - 1)
    free(abuffer - 1)
    Free_Work_Data(work)
  Close_DBX(dbx1)
  if ISTWO: Close_DBX(dbx2)
  exit(0)
