# for C++ define  CC = g++
CC = g++
CFLAGS = -g -Wall -fPIC -m64
LFLAGS = -L$(ROOTSYS)/lib -L$(ALICE_ROOT)/lib
INC =	-I$(ROOTSYS)/include  -I$(ALICE_ROOT)/include -I./
TGT =	libAlg.so
DICT=	AlgDict.cxx
DICTO=	AlgDict.o

SRC = 	AliAlgAux.cxx AliAlgDetTPC.cxx AliAlgSens.cxx AliAlgSensTRD.cxx  \
	AliAlgDet.cxx AliAlgDetTRD.cxx AliAlgSensITS.cxx \
	AliAlgSteer.cxx AliAlgDetITS.cxx  AliAlgMPRecord.cxx  AliAlgSensTOF.cxx \
	AliAlgTrack.cxx AliAlgDetTOF.cxx  AliAlgPoint.cxx AliAlgSensTPC.cxx \
	AliAlgVol.cxx AliAlgVtx.cxx AliAlgRes.cxx \
	Mille.cxx


HDR =	$(SRC:.cxx=.h)

OBJ = 	$(SRC:.cxx=.o)


.PHONY: depend clean

all: 	$(TGT)
	@echo creating libAlg.so

$(TGT):	$(OBJ) $(DICTO)
	$(CC) $(CFLAGS)  -shared -o $(TGT) $(OBJ) $(DICTO) `root-config --ldflags` $(LFLAGS)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

%.o : %.cxx
	$(CC) $(CFLAGS) $(INC) -c $<  -o $@
	$(CC) -MM $(CFLAGS) $(INC) -c $*.cxx > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm *.o *~ *.so *.d AlgDict.{h,cxx}

$(DICT): $(HDR) AlgLinkDef.h
	rootcint -f $@ -c $(INC) $(HDR) $^


depend: $(SRC)
	makedepend $(INC) $^

# DO NOT DELETE THIS LINE -- make depend needs it
