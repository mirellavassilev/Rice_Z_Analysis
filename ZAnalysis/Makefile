ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

MKDIR_BIN=mkdir -p $(PWD)/bin


all: mkdirBin bin/Z_EE_Channel.exe

mkdirBin:
	$(MKDIR_BIN)


bin/Z_EE_Channel.exe: src/Z_EE_Channel.C
	$(CXX) $(CXXFLAGS) $(ROOT) -o bin/Z_EE_Channel.exe src/Z_EE_Channel.C

clean:
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bin/*.exe
	rmdir bin
