CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET = main
SRCS = src/pdb_cpp/main.cpp src/pdb_cpp/Model.cpp  src/pdb_cpp/Coor.cpp src/pdb_cpp/format/pdb.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) $(OBJS)