CC=gcc
GROMACS=/usr/local/gromacs
VGRO=5
INCLUDE=include
SRC=src
BUILD=build
INSTALL=/usr/local/bin

ifeq ($(VGRO),5)
INCGRO=-I$(GROMACS)/include/ -I$(GROMACS)/include/gromacs/utility -I$(GROMACS)/include/gromacs/fileio -I$(GROMACS)/include/gromacs/commandline -I$(GROMACS)/include/gromacs/legacyheaders
LINKGRO=-L$(GROMACS)/lib/i386-linux-gnu
LIBGRO=-lgromacs
DEFV5=-D GRO_V5
else
INCGRO=-I$(GROMACS)/include/gromacs
LINKGRO=-L$(GROMACS)/lib
LIBGRO=-lgmx
endif

.PHONY: install clean

$(BUILD)/tr_tessellator: $(BUILD)/tr_tessellator.o $(BUILD)/trt_tessellation.o
	$(CC) $(CFLAGS) -o $(BUILD)/tr_tessellator $(BUILD)/tr_tessellator.o $(BUILD)/trt_tessellation.o $(LINKGRO) $(LIBGRO)

install: $(BUILD)/tr_tessellator
	install $(BUILD)/tr_tessellator $(INSTALL)

$(BUILD)/tr_tessellator.o: $(SRC)/tr_tessellator.c $(INCLUDE)/trt_tessellation.h
	$(CC) $(CFLAGS) -o $(BUILD)/tr_tessellator.o -c $(SRC)/tr_tessellator.c $(DEFV5) -I$(INCLUDE) $(INCGRO)

$(BUILD)/trt_tessellation.o: $(SRC)/trt_tessellation.c $(INCLUDE)/trt_tessellation.h
	$(CC) $(CFLAGS) -o $(BUILD)/trt_tessellation.o -c $(SRC)/trt_tessellation.c $(DEFV5) -I$(INCLUDE) $(INCGRO)

clean:
	rm -f $(BUILD)/tr_tessellator.o $(BUILD)/tr_tessellator