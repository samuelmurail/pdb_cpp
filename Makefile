CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET = main
SRCS = src/pdb_cpp/select.cpp src/pdb_cpp/Model.cpp  src/pdb_cpp/Coor.cpp src/pdb_cpp/sequence.cpp src/pdb_cpp/sequence_align.cpp src/pdb_cpp/format/pdb.cpp src/pdb_cpp/data/residue.cpp src/pdb_cpp/main.cpp 
OBJS = $(SRCS:.cpp=.o)

scratch: clean $(TARGET)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) $(OBJS)