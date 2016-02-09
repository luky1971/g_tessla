CC ?= gcc
CFLAGS += -std=c99 -O3 -DLLT_BENCH
# CFLAGS += -std=c99 -DLLT_DEBUG -g
# CFLAGS += -std=c99 -g

GROMACS = /usr/local/gromacs
VGRO = 5

GKUT = extern/gkut
PRED = extern/predicates

INCLUDE = include
SRC = src
BUILD = build
INSTALL = /usr/local/bin

LIBS = -lm

ifeq ($(VGRO),5)
INCGRO = -I$(GROMACS)/include/ \
	-I$(GROMACS)/include/gromacs/utility \
	-I$(GROMACS)/include/gromacs/fileio \
	-I$(GROMACS)/include/gromacs/commandline \
	-I$(GROMACS)/include/gromacs/legacyheaders
LINKGRO = -L$(GROMACS)/lib/x86_64-linux-gnu
LIBGRO = -lgromacs
DEFV5 = -D GRO_V5
else
INCGRO = -I$(GROMACS)/include/gromacs
LINKGRO = -L$(GROMACS)/lib
LIBGRO = -lgmx
endif

ifneq ($(PARALLEL),0)
CFLAGS += -fopenmp
endif

MCFLAGS ='
MCFLAGS +=$(CFLAGS)
MCFLAGS +='

.PHONY: install clean

$(BUILD)/lltessellator: $(BUILD)/lltessellator.o $(BUILD)/llt_tri.o $(BUILD)/llt_grid.o $(BUILD)/delaunay_tri.o
	make CC=$(CC) CFLAGS=$(MCFLAGS) GROMACS=$(GROMACS) VGRO=$(VGRO) -C $(GKUT) \
	&& make CC=$(CC) -C $(PRED) \
	&& $(CC) $(CFLAGS) -o $(BUILD)/lltessellator $(BUILD)/lltessellator.o $(BUILD)/llt_tri.o $(BUILD)/llt_grid.o $(BUILD)/delaunay_tri.o \
	$(GKUT)/build/gkut_io.o $(GKUT)/build/gkut_log.o $(PRED)/predicates.o $(LINKGRO) $(LIBGRO) $(LIBS)

install: $(BUILD)/lltessellator
	install $(BUILD)/lltessellator $(INSTALL)

$(BUILD)/lltessellator.o: $(SRC)/lltessellator.c $(INCLUDE)/llt_grid.h $(INCLUDE)/llt_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/lltessellator.o -c $(SRC)/lltessellator.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include

$(BUILD)/llt_tri.o: $(SRC)/llt_tri.c $(INCLUDE)/llt_tri.h $(INCLUDE)/delaunay_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/llt_tri.o -c $(SRC)/llt_tri.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(PRED)

$(BUILD)/llt_grid.o: $(SRC)/llt_grid.c $(INCLUDE)/llt_grid.h
	$(CC) $(CFLAGS) -o $(BUILD)/llt_grid.o -c $(SRC)/llt_grid.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include

$(BUILD)/delaunay_tri.o: $(SRC)/delaunay_tri.c $(INCLUDE)/delaunay_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/delaunay_tri.o -c $(SRC)/delaunay_tri.c -I$(INCLUDE) -I$(PRED)

clean:
	make clean -C $(GKUT) \
	&& make clean -C $(PRED) \
	&& rm -f $(BUILD)/*.o $(BUILD)/lltessellator
