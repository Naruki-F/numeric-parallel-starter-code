INCLUDE_DIRS = 
LIB_DIRS = 
CC = gcc

CDEFS=
CFLAGS= -O3 -Wall $(INCLUDE_DIRS) $(CDEFS)
LIBS=-lpthread

PRODUCT=sharpen_grid sharpen

CFILES= sharpen_grid.c sharpen.c


SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~ *sharp*.ppm
	-rm -f ${PRODUCT}

sharpen_grid:	sharpen_grid.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ sharpen_grid.o $(LIBS)

sharpen:	sharpen.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ sharpen.o $(LIBS)

depend:

.c.o:
	$(CC) $(CFLAGS) -c $<