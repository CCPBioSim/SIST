#######################
# Makefile for qsidd  #
#######################

OBJECTS = \
	G_x.o\
	SIDD_1R.o\
	SIDD_2R.o\
	SIDD_3R.o\
	SIDD_4R.o\
	SIDD_Base.o\
	stat_1R.o\


APP = qsidd
SRC = qsidd.c
OBJ = qsidd.o


DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ) $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(APP) $(CFLAGS) $(OBJ) $(OBJECTS) -lm

clean:
	rm -f *.o $(APP)

depend: $(OBJECTS:.o=.cpp)
	$(CXX) -MM $^ > $@

test:
	./qsidd example.fasta 

tar:
	rm -rf /tmp/$(APP)
	mkdir /tmp/$(APP)
	cp -r * /tmp/$(APP)
	cd /tmp/$(APP); make clean
	cd /tmp; tar -zcvf $(APP)-$(DATE).tar.gz $(APP)
	rm -rf /tmp/$(APP)


#################
# Architectures #
#################

gcc:
	make $(APP) CC="$(CXX)" CPPFLAGS="-O2 -Wall"


###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CPPFLAGS) -c -o $@ $<

################
# Dependancies #
################

include depend

