# ============================================================================
#  Makefile HBTeg                             Chris Plumberg, October 31, 2018
# ============================================================================
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

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)
INCDIR		=	include
SRCDIR		=	src

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	shifter.e
endif

SRC			=	main.cpp \
				$(SRCDIR)/shifter.cpp \
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

OBJDIR			=	.
OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
					$(basename $(SRC))))
#OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
#					$(notdir $(basename $(SRC)))))
TARGET			=	$(MAIN)
INSTPATH		=	..

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install target

all:		mkobjdir $(TARGET)

target:		mkobjdir $(TARGET)

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

