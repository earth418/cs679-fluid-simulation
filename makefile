LDFLAGS_COMMON = -lGL -lGLU -lglut -lstdc++
CFLAGS_COMMON = -c -Wall -I./ -O3
COMPILER = g++

# calls:
CC         = g++
CFLAGS     = ${CFLAGS_COMMON} -O3
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = HW4

SOURCES    = HW4.cpp \
	     FIELD_2D.cpp \
	     FLUID_2D.cpp \
	     FLUID_2D_BOUNDED.cpp \
	     FLUID_2D_PERIODIC.cpp

OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
