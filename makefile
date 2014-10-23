EXTERNAL_OBJS = triangle.o
SRCS = Tricall_Pre_Mesh_Gen.cc Tricall_Utility.cc Data_Transfer_Object.cc
INCLUDES = Tricall_Pre_Mesh_Gen.h Tricall_Utility.h Data_Transfer_Object.h
OBJS = ${SRCS:.cc=.o}
TARG = pdt 
CC = mpicxx 
BOOST_ROOT=/usr/local/boost-1.53.0
STAPL_ROOT=/home/akansha/stapl
CCFLAGS = -O0 -g -std=c++11 -D_STAPL -DNO_STAPL_MAIN\
	  -I$(STAPL_ROOT)/tools/libstdc++/4.7.2 \
	  -I$(STAPL_ROOT)/tools \
	  -I$(STAPL_ROOT) \
	  -I$(BOOST_ROOT)/include
LIB = -L$(STAPL_ROOT)/lib  -lstapl -lstapl_rt -L$(BOOST_ROOT)/lib -lboost_serialization -lboost_system

.cc.o:
	$(CC) $(CCFLAGS) -c $<

default:
	@$(MAKE) CCFLAGS="$(CCFLAGS) $(OPT_FLAGS)" TARG="$(TARG)" $(TARG)

clean:
	/bin/rm -rf $(OBJS) $(TARG) core

$(TARG): $(INCLUDES) $(OBJS)
	$(CC) -o $(TARG) $(CCFLAGS) $(OBJS) $(EXTERNAL_OBJS) $(LIB)
