CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET = main
SRCS = src/pdb_cpp/_core/select.cpp src/pdb_cpp/_core/Model.cpp src/pdb_cpp/_core/Coor.cpp src/pdb_cpp/_core/sequence.cpp src/pdb_cpp/_core/align.cpp src/pdb_cpp/_core/seq_align.cpp src/pdb_cpp/_core/format/pdb.cpp src/pdb_cpp/_core/data/residue.cpp src/pdb_cpp/_core/main.cpp
OBJS = $(SRCS:.cpp=.o)

scratch: clean $(TARGET)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) $(OBJS)