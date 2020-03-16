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

#OBJECTS				=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
#						$(basename $(SRC))))
#OBJECTS					=	$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,\
#						$(sort $(wildcard $(SRCDIR)/*.cpp)))
OBJECTS				=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
						$(notdir $(basename $(SRC)))))
OBJECTS				+=	$(addsuffix $O, $(basename $(MAINSRC)))
#OBJECTS				=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
#						$(notdir $(basename $(SRC)))))
TARGET				=	$(MAIN)
TARGET_LIBRARY		=	$(LIBDIR)/libshifter.a
TARGET_DYN_LIBRARY	=	$(LIBDIR)/libshifter.so
INSTPATH			=	..

# --------------- Pattern rules -------------------

%.o: %.cpp
	$(CC) $(CFLAGS) -c -fPIC $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c -fPIC $< -o $@

$(SRCDIR)/%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

$(TARGET):	$(OBJECTS)	
	-@mkdir -p $(LIBDIR)
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 

$(TARGET_LIBRARY): $(OBJECTS)
	-@mkdir -p $(LIBDIR)
	ar cru $@ $^
	ln -sf `readlink -e $@` ~/$(TARGET_LIBRARY)

$(TARGET_DYN_LIBRARY): $(OBJECTS)
	-@mkdir -p $(LIBDIR)
	$(CC) -shared $^ -o $@
	ln -sf `readlink -e $@` ~/$(TARGET_DYN_LIBRARY)

# -------------------------------------------------

.PHONY:		all mkdirs clean distclean install target lib

all:		mkdirs $(TARGET) $(TARGET_LIBRARY) $(TARGET_DYN_LIBRARY)

target:		mkdirs $(TARGET)

lib:		mkdirs $(TARGET_LIBRARY) $(TARGET_DYN_LIBRARY)

help:
	@grep '^##' GNUmakefile

mkdirs:	
	-@mkdir -p $(OBJDIR)
	-@mkdir -p $(LIBDIR)

clean:		
	-rm -f $(OBJECTS)

distclean:	
	-rm -f $(TARGET)
	-rm -f $(TARGET_LIBRARY)
	-rm -f $(TARGET_DYN_LIBRARY)
	-rm -f $(OBJECTS)

install:	$(TARGET)
	cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/Arsenal.cpp: $(INCDIR)/Arsenal.h
$(SRCDIR)/ParameterReader.cpp: $(INCDIR)/ParameterReader.h $(INCDIR)/Arsenal.h
$(SRCDIR)/shifter.cpp: $(INCDIR)/shifter.h $(INCDIR)/ParameterReader.h \
						$(INCDIR)/Arsenal.h $(INCDIR)/ParticleRecord.h
./main.cpp: $(INCDIR)/shifter.h $(INCDIR)/ParameterReader.h $(INCDIR)/ParticleRecord.h

