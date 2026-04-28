CC ?=			gcc
CXX ?=		g++
CFLAGS ?=		-g -Wall -O3
CXXFLAGS ?=	$(CFLAGS)
CPPFLAGS=
PROG ?=		minipileup2
OBJS ?=		bedidx.o

# htslib: set HTSLIB_DIR to override (e.g. make HTSLIB_DIR=/path/to/htslib)
# otherwise pkg-config is tried, then system default
ifdef HTSLIB_DIR
  INCLUDES=	-I$(HTSLIB_DIR)
  HTSLIB_LIBS=	-L$(HTSLIB_DIR) -lhts -Wl,-rpath,$(HTSLIB_DIR)
else
  INCLUDES=	$(shell pkg-config --cflags htslib 2>/dev/null)
  HTSLIB_LIBS=	$(shell pkg-config --libs htslib 2>/dev/null || echo "-lhts")
endif

LIBS ?=		$(HTSLIB_LIBS) -lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

minipileup2:$(OBJS) pileup.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

bedidx.o: ksort.h kseq.h khash.h
