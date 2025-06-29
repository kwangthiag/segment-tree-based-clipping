# How to make with debugging on
# make clean; make D=0 S=0
# How to make withoout debugging
# make clean; make

# CC = g++
# #CFLAGS = -Wall -c -std=c++0x
# CFLAGS = -c -O2 -std=c++0x
# DEBUG = -g

CC = pgc++
# #CFLAGS = -Wall -c -std=c++0x
# CFLAGS = -c -acc
# DEBUG =  -fast -fastsse -Mipa=fast -g  
# GOPT = 
# OPT = 
# CFLAGS = -c -acc
LIBCUDA = -L$(HOME)/nvhpc/Linux_x86_64/25.5/cuda/12.9/lib64


# DEBUG =  -Minfo=all
DEBUG = -fast -fastsse -Mipa=fast -g
GOPT = -Xptxas -O2
OPT = -O2
CFLAGS = -c -acc -mp -O2 -Ilib/PAM/include -Ilib/parlaylib/include -Ip_seg_tree/segment_tree

DEBUG += -DTIME # TIME flag is used to extract run times in the functions
DEBUG += -DTIME_L2 # TIME_L2 flag is used to extract run times inside parallel functions
ifdef D
    DEBUG += -DDG -DDT_INFO -DDG2 -DDG3
endif

ifdef S
    DEBUG += -DSAVE
endif

ifdef gpu 
	DEBUG += -DGPU
else
	DEBUG += -DCPU
endif

# default compilation
all: cpu  

cpu: util.o thrust_func.o oacc_segment_tree.o  clip.o 
	# $(CC) -acc -O2 -o bin/seg bin/util.o thrust_func.o oacc_segment_tree.o bin/clip.o 
	$(CC) -acc $(OPT) -mp -o bin/seg bin/util.o bin/oacc_segment_tree.o bin/clip.o $(LIBCUDA) -lcudart

gpu: util.o thrust_func.o oacc_segment_tree_gpu.o  clip_gpu.o 
	$(CC) -acc $(OPT) -o bin/seg bin/util.o bin/thrust_func.o bin/oacc_segment_tree_gpu.o bin/clip_gpu.o  $(LIBCUDA) -lcudart
	
pam: util.o thrust_func.o pam_segment_tree.o clip.o
	$(CC) -acc $(OPT) -mp -o bin/seg bin/util.o bin/pacc_segment_tree.o bin/clip.o $(LIBCUDA) -lcudart

clip.o: src/clip.cpp
	# $(CC) $(CFLAGS) $(DEBUG) -o clip.o -c clip.cpp 
	$(CC) $(CFLAGS) $(DEBUG) -o bin/clip.o -c src/clip.cpp -fopenmp
	# $(CC) $(CFLAGS) $(DEBUG) -o bin/clip.o -c srcclip.cpp

tree: driver.o util.o segment_tree.o
	$(CC) $(OPT)  -o seg bin/util.o segment_tree.o driver.o
	
driver.o:  src/driver.cpp
	$(CC) $(CFLAGS) $(DEBUG) src/driver.cpp
		
util.o: p_seg_tree/util.cpp
	$(CC) $(CFLAGS) $(DEBUG) -o bin/util.o p_seg_tree/util.cpp
	
segment_tree.o: p_seg_tree/segment_tree.cpp
	$(CC) $(CFLAGS) $(DEBUG) p_seg_tree/segment_tree.cpp 

oacc_segment_tree.o: p_seg_tree/oacc_segment_tree.cpp
	# $(CC) -g $(CFLAGS) $(DEBUG) p_seg_tree/oacc_segment_tree.cpp -fopenmp	 
	$(CC) $(CFLAGS) $(DEBUG) -o bin/oacc_segment_tree.o p_seg_tree/oacc_segment_tree.cpp -fopenmp	 
	# $(CC)  $(CFLAGS) $(DEBUG) p_seg_tree/oacc_segment_tree.cpp	


oacc_segment_tree_gpu.o: p_seg_tree/oacc_segment_tree.cpp
	# $(CC) -g $(CFLAGS) $(DEBUG) p_seg_tree/oacc_segment_tree.cpp -fopenmp	 
	$(CC) $(CFLAGS) $(DEBUG) -o bin/oacc_segment_tree_gpu.o p_seg_tree/oacc_segment_tree.cpp -fopenmp	

clip_gpu.o: src/clip.cpp
	# $(CC) $(CFLAGS) $(DEBUG) -o bin/clip.o -c clip.cpp 
	$(CC) $(CFLAGS) $(DEBUG) -o bin/clip_gpu.o -c src/clip.cpp -fopenmp

thrust_func.o: p_seg_tree/thrust_func.cu
	nvcc $(GOPT) -o bin/thrust_func.o -c p_seg_tree/thrust_func.cu -ccbin /usr/bin/gcc-12 -allow-unsupported-compiler 

# pam_segment_tree.o: p_seg_tree/pam_segment_tree.cpp
#     $(CC) $(CFLAGS) $(DEBUG) -o bin/pam_segment_tree.o p_seg_tree/pam_segment_tree.cpp -fopenmp   
pam_segment_tree.o: p_seg_tree/pam_segment_tree.cpp p_seg_tree/segment_tree/pam_segment_tree.h
	$(CC) $(CFLAGS) $(DEBUG) -c p_seg_tree/pam_segment_tree.cpp -o $@

# polyclip_time.o: lib/optimizedFostersAlgorithm/polyclip_time.cpp
# 	$(CC) -g $(CFLAGS) $(DEBUG) lib/optimizedFostersAlgorithm/polyclip_time.cpp

clean:
	rm bin/*.o	
	rm bin/seg
