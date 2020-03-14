# =============================================================
#  Makefile                             Christopher J. Plumberg
# =============================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3 -std=c++11 -g

RM				=	rm -f
O               =	.o
LDFLAGS         =	$(CFLAGS)
SYSTEMFILES     =	$(SRCGNU)
INCDIR			=	include
SRCDIR			=	src
LIBDIR			=	lib
OBJDIR			=	obj

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	shifter.e
endif

MAINSRC		=	main.cpp

SRC			=	$(SRCDIR)/shifter.cpp \
				$(SRCDIR)/ParameterReader.cpp \
				$(SRCDIR)/Arsenal.cpp \
				$(SRCDIR)/ParticleRecord.cpp \
				$(SRCDIR)/FourVector.cpp

INC			= 	$(INCDIR)/random_events.h \
				$(INCDIR)/shifter.h \
				$(INCDIR)/ParameterReader.h \
				$(INCDIR)/Arsenal.h \
				$(INCDIR)/ParticleRecord.h \
				$(INCDIR)/FourVector.h

# -------------------------------------------------

#OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
#					$(basename $(SRC))))
OBJECTS				=	$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,\
					$(sort $(wildcard $(SRCDIR)/*.cpp)))
#OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
#					$(notdir $(basename $(SRC)))))
TARGET			=	$(MAIN)
TARGET_LIBRARY	=	$(LIBDIR)/libshifter.a
INSTPATH		=	..

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/main.o: $(MAINSRC)
	$(CC) $(CFLAGS) -c $< -o $@
	OBJECTS+=$(OBJDIR)/main.o

$(LIBDIR)/libshifter.a: $(OBJECTS)
	#rm -f $(LIBDIR)/libpythia8$(LIB_SUFFIX)
	ar cru $@ $^

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install target lib

all:		mkobjdir $(TARGET) $(TARGET_LIBRARY)

target:		mkobjdir $(TARGET)

lib:		mkobjdir $(TARGET_LIBRARY)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 

clean:		
		-rm -f $(OBJECTS)

distclean:	
		-rm -f $(TARGET)
		-rm -f $(OBJECTS)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/Arsenal.cpp: $(INCDIR)/Arsenal.h
$(SRCDIR)/ParameterReader.cpp: $(INCDIR)/ParameterReader.h $(INCDIR)/Arsenal.h
$(SRCDIR)/shifter.cpp: $(INCDIR)/shifter.h $(INCDIR)/ParameterReader.h \
						$(INCDIR)/Arsenal.h $(INCDIR)/ParticleRecord.h
./main.cpp: $(INCDIR)/shifter.h $(INCDIR)/ParameterReader.h $(INCDIR)/ParticleRecord.h

