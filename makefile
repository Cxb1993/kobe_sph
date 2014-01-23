FC = ifort
FFLAGS = -O0
#FFLAGS = -O3 -opt-report

SRC = particle.f90 main.f90

TARGET = sph

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f $(TARGET) *.o *.mod
