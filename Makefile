#
#
#   General Makefile to the POLYS programme (UNIX using X11)
#
#
#   System dependent names, flags and options:
#   ==========================================
#
#   CC=		cc       most UNIX systems
#   CC=		xlc      IBM AIX systems
#   CC=		gcc      GNU C compiler
#
#   CFLAGS=	-O       most UNIX systems
#               -Q       creates inline information on the xlc compiler
#
#   PXLIBS=	-lX11 -lm          -these 2 flags are always required.
#   PXLIBS=	-lX11 -lm -lbsd    -required on many systems.
#                                   e.g.: Silicon Graphics IRIX, IBM AIX
#

CC=		gcc
CFLAGS=		-O
# CFLAGS=		-g
# CFLAGS=		-O -Q
# CFLAGS=		-g -fullwarn -prototypes

PXLIBS=		-L/usr/X11R6/lib/ -lX11 -lm -lbsd
PMLIBS=		-lm -ll
POLYSDIR=	~/polys

POLYS=		polys
PHIPSI=         phipsi
UNIMERGE=       unimerge
MAPAREA=        maparea
MERGEOR=        mergeor
BOLTZ=          boltz
TG=             trans.gauche
TRANS=          translate

IO_OBJ=         sinput.o \
                fileoper.o

MM_OBJ=         lex.yy.o \
                polypars.tab.o \
                action.o \
                bondtab.o \
                brent.o \
                calculat.o \
                coord.o \
                driver.o \
                doubleh.o \
                expand.o \
                f1dim.o \
                frac2carte.o \
                generate.o \
                helix.o \
                initiate.o \
                intcoord.o \
                invmat.o \
                konjug.o \
                glyclink.o \
                linmin.o \
                linsch.o \
                matrix.o \
                minimize.o \
                mnbrak.o \
                orient.o \
                palloc.o \
                phelp.o \
                phimin.o \
                polys.o \
                potener.o \
                powell.o \
                qsort.o \
                r_charmm.o \
                r_csd.o \
                r_dat.o \
                r_mm3.o \
                r_mono.o \
                r_pdb.o \
                r_polys.o \
                r_sybyl.o \
                random.o \
                readcoor.o \
                readpara.o \
                rotate.o \
                timer.o \
                vbond.o \
                vtheta.o \
                vphi.o \
                vcon.o \
                vnonbond.o \
                writcoor.o 


POLYS_OBJ=	${IO_OBJ} ${MM_OBJ}

PHIPSI_OBJ=     ${IO_OBJ} phipsi.o

UNIMERGE_OBJ=	${IO_OBJ} unimerge.o

MAPAREA_OBJ=	${IO_OBJ} maparea.o

MERGEOR_OBJ=	${IO_OBJ} mergeor.o

BOLTZ_OBJ=	${IO_OBJ} boltz.o

TG_OBJ=		${IO_OBJ} trans.gauche.o

TRANS_OBJ=	${IO_OBJ} translate.o

all:		polys phipsi unimerge maparea mergeor boltz trans.gauche translate

lex.yy.o:       polypars.tab.h  defs.h typs.h polys.h
y.tab.o:                 defs.h typs.h polys.h
brent.o:
driver.o:       polys.h  extern.h
doubleh.o:      defs.h typs.h extern.h polys.h

f1dim.o:        defs.h
generate.o:     polys.h  extern.h
intcoord.o:     polys.h  extern.h matrix.h
konjug.o:       polys.h  extern.h
helix.o:        defs.h typs.h extern.h polys.h
initiate.o:     polys.h  extern.h
glyclink.o:     polys.h  extern.h
linmin.o:
linsch.o:       defs.h
matrix.o:       defs.h typs.h polys.h
minimize.o:     polys.h  extern.h
mnbrak.o: 
orient.o:       polys.h  extern.h
palloc.o:      
phimin.o:       polys.h  extern.h
polys.o:        defs.h typs.h vars.h   polys.h
qsort.o:        defs.h
random.o:      
potener.o:      defs.h typs.h extern.h matrix.h 
powell.o:       polys.h defs.h typs.h extern.h
readcoor.o:     defs.h typs.h extern.h fileoper.h
readpara.o:     defs.h typs.h extern.h polys.h
rotate.o:       defs.h typs.h extern.h polys.h
timer.o:	
vbond.o:        defs.h typs.h extern.h sinput.h fileoper.h matrix.h
vtheta.o:       defs.h typs.h extern.h sinput.h fileoper.h matrix.h
vphi.o:         defs.h typs.h extern.h sinput.h fileoper.h matrix.h
vcon.o:         defs.h typs.h extern.h sinput.h fileoper.h matrix.h
vnonbond.o:     defs.h typs.h extern.h          fileoper.h
writcoor.o:     defs.h typs.h extern.h          fileoper.h

polys:		${POLYS_OBJ}
		${CC} ${POLYS_OBJ} ${PMLIBS} -o ${POLYS}

phipsi:		${PHIPSI_OBJ}
		${CC} ${PHIPSI_OBJ} ${PMLIBS} -o ${PHIPSI}

unimerge:	${UNIMERGE_OBJ}
		${CC} ${UNIMERGE_OBJ} ${PMLIBS} -o ${UNIMERGE}

maparea:	${MAPAREA_OBJ}
		${CC} ${MAPAREA_OBJ} ${PMLIBS} -o ${MAPAREA}

mergeor:	${MERGEOR_OBJ}
		${CC} ${MERGEOR_OBJ} ${PMLIBS} -o ${MERGEOR}

boltz:    	${BOLTZ_OBJ}
		${CC} ${BOLTZ_OBJ} ${PMLIBS} -o ${BOLTZ}

trans.gauche:  	${TG_OBJ}
		${CC} ${TG_OBJ} ${PMLIBS} -o ${TG}

translate:  	${TRANS_OBJ}
		${CC} ${TRANS_OBJ} ${PMLIBS} -o ${TRANS}

clean:		
		rm -f *.o
