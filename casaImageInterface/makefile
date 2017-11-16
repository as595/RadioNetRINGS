#------------------------------------------------------------------------------
# Update this constant to the casacore folder path.
CASACORE=/home/cjs105/Documents/code/casacore-master

# Nothing after this point should have to be altered.
MYPROGRAM=cppTestCasaImageInterface
CODEFOLDER=./
CASAIMAGEINTERFACE=cppCasaImageInterface
CASACOREBUILD=$(CASACORE)/build

# folders for the casa run-time libraries.
CASA_DIR=$(CASACOREBUILD)/casa
LATTICES_DIR = $(CASACOREBUILD)/lattices
COORDINATES_DIR = $(CASACOREBUILD)/coordinates
IMAGES_DIR = $(CASACOREBUILD)/images

# the libraries themselves.
CASA=$(CASA_DIR)/libcasa_casa.so #.$(CASACORE_VERSION)
LATTICES=$(LATTICES_DIR)/libcasa_lattices.so #.$(CASACORE_VERSION)
COORDINATES=$(COORDINATES_DIR)/libcasa_coordinates.so #.$(CASACORE_VERSION)
IMAGES=$(IMAGES_DIR)/libcasa_images.so #.$(CASACORE_VERSION)
CC=g++

#------------------------------------------------------------------------------

all: $(MYPROGRAM)

$(MYPROGRAM): $(CODEFOLDER)/$(MYPROGRAM).cpp $(CASAIMAGEINTERFACE).o $(CASA) $(LATTICES) $(COORDINATES) $(IMAGES)
	$(CC) $^ -o $@ -I $(CASACORE) -I $(CASACOREBUILD) -Wl,-rpath,$(CASA_DIR):$(LATTICES_DIR):$(COORDINATES_DIR):$(IMAGES_DIR)

$(CASAIMAGEINTERFACE).o: $(CODEFOLDER)/$(CASAIMAGEINTERFACE).cpp
	$(CC) -c $^ -o $@ -I $(CASACORE) -I $(CASACOREBUILD)

clean:
	rm -f $(MYPROGRAM)
