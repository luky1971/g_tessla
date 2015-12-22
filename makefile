CC ?= gcc
CFLAGS += -O3

GROMACS = /usr/local/gromacs
VGRO = 5

GKUT = extern/gkut
TRI = extern/triangle

INCLUDE = include
SRC = src
BUILD = build
INSTALL = /usr/local/bin

ifeq ($(VGRO),5)
INCGRO = -I$(GROMACS)/include/ \
	-I$(GROMACS)/include/gromacs/utility \
	-I$(GROMACS)/include/gromacs/fileio \
	-I$(GROMACS)/include/gromacs/commandline \
	-I$(GROMACS)/include/gromacs/legacyheaders
LINKGRO = -L$(GROMACS)/lib/i386-linux-gnu
LIBGRO = -lgromacs
DEFV5 = -D GRO_V5
else
INCGRO = -I$(GROMACS)/include/gromacs
LINKGRO = -L$(GROMACS)/lib
LIBGRO = -lgmx
endif

.PHONY: install clean

$(BUILD)/lltessellator: $(BUILD)/lltessellator.o $(BUILD)/llt_tri.o $(BUILD)/llt_grid.o
	make GROMACS=$(GROMACS) VGRO=$(VGRO) -C $(GKUT) \
	&& make trilibrary -C $(TRI) \
	&& $(CC) $(CFLAGS) -o $(BUILD)/lltessellator $(BUILD)/lltessellator.o $(BUILD)/llt_tri.o $(BUILD)/llt_grid.o \
	$(GKUT)/build/gkut_io.o $(GKUT)/build/gkut_log.o $(TRI)/triangle.o $(LINKGRO) $(LIBGRO)

install: $(BUILD)/lltessellator
	install $(BUILD)/lltessellator $(INSTALL)

$(BUILD)/lltessellator.o: $(SRC)/lltessellator.c $(INCLUDE)/llt_grid.h $(INCLUDE)/llt_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/lltessellator.o -c $(SRC)/lltessellator.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(TRI)

$(BUILD)/llt_tri.o: $(SRC)/llt_tri.c $(INCLUDE)/llt_tri.h
	$(CC) $(CFLAGS) -o $(BUILD)/llt_tri.o -c $(SRC)/llt_tri.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(TRI)

$(BUILD)/llt_grid.o: $(SRC)/llt_grid.c $(INCLUDE)/llt_grid.h
	$(CC) $(CFLAGS) -o $(BUILD)/llt_grid.o -c $(SRC)/llt_grid.c \
	$(DEFV5) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include

clean:
	make clean -C $(GKUT) && make distclean -C $(TRI) && rm -f $(BUILD)/*.o $(BUILD)/lltessellator