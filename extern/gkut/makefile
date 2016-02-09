CC ?= gcc
CFLAGS += -O3
GROMACS = /usr/local/gromacs
VGRO = 5
INCLUDE = include
SRC = src
BUILD = build

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

.PHONY: all

all: $(BUILD)/gkut_io.o $(BUILD)/gkut_log.o

$(BUILD)/gkut_io.o: $(SRC)/gkut_io.c $(INCLUDE)/gkut_io.h
	$(CC) $(CFLAGS) -o $(BUILD)/gkut_io.o -c $(SRC)/gkut_io.c $(DEFV5) -I$(INCLUDE) $(INCGRO)

$(BUILD)/gkut_log.o: $(SRC)/gkut_log.c $(INCLUDE)/gkut_log.h
	$(CC) $(CFLAGS) -o $(BUILD)/gkut_log.o -c $(SRC)/gkut_log.c $(DEFV5) -I$(INCLUDE) $(INCGRO)

clean:
	rm -f $(BUILD)/*.o
