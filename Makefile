FC = mpif77
OBJSUF = .o

objects.for := $(wildcard *.for)
objects.for.obj:= $(objects.for:.for=$(OBJSUF))

FCFLAGS =
# -pg

all: moto.exe
moto.exe: $(objects.for.obj)
	$(FC) $(FCFLAGS) -o $@ $^

$(objects.for.obj):%$(OBJSUF): %.for
	$(FC) -c $^

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o

veryclean: clean
	rm -f *~ $(PROGRAMS)
