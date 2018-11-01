include $(FSLCONFDIR)/default.mk
CXX=g++
PROJNAME = filtershift

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

OBJS = filtershift.o Window.o bessel.o

LIBS = -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm -lutils -lz 

XFILES = filtershift

all:	${XFILES}

filtershift:${OBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${OBJS} ${LIBS}

