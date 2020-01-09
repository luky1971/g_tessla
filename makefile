CC = gcc
CFLAGS += -std=c99 -O3
# CFLAGS += -std=c99 -O3 -DGTA_BENCH
# CFLAGS += -std=c99 -g -DGTA_DEBUG
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

$(BUILD)/g_tessla: $(BUILD)/g_tessla.o $(BUILD)/gta_tri.o $(BUILD)/gta_grid.o $(BUILD)/delaunay_tri.o
	make CC=$(CC) CFLAGS=$(MCFLAGS) GROMACS=$(GROMACS) VGRO=$(VGRO) -C $(GKUT) \
	&& make CC=$(CC) -C $(PRED) \
	&& $(CC) $(CFLAGS) -o $(BUILD)/g_tessla $(BUILD)/g_tessla.o $(BUILD)/gta_tri.o $(BUILD)/gta_grid.o $(BUILD)/delaunay_tri.o \
	$(GKUT)/build/gkut_io.o $(GKUT)/build/gkut_log.o $(PRED)/predicates.o $(LINKGRO) $(LIBGRO) $(LIBS)

install: $(BUILD)/g_tessla
	install $(BUILD)/g_tessla $(INSTALL)

$(BUILD)/g_tessla.o: $(SRC)/g_tessla.c $(INCLUDE)/gta_grid.h $(INCLUDE)/gta_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/g_tessla.o -c $(SRC)/g_tessla.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include

$(BUILD)/gta_tri.o: $(SRC)/gta_tri.c $(INCLUDE)/gta_tri.h $(INCLUDE)/delaunay_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/gta_tri.o -c $(SRC)/gta_tri.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(PRED)

$(BUILD)/gta_grid.o: $(SRC)/gta_grid.c $(INCLUDE)/gta_grid.h
	$(CC) $(CFLAGS) -o $(BUILD)/gta_grid.o -c $(SRC)/gta_grid.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include

$(BUILD)/delaunay_tri.o: $(SRC)/delaunay_tri.c $(INCLUDE)/delaunay_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/delaunay_tri.o -c $(SRC)/delaunay_tri.c -I$(INCLUDE) -I$(PRED)

clean:
	make clean -C $(GKUT) \
	&& make clean -C $(PRED) \
	&& rm -f $(BUILD)/*.o $(BUILD)/g_tessla
