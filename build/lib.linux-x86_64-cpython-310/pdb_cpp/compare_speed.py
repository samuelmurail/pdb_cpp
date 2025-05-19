import pdb_cpp
import time
import pdb_numpy


start_time = time.time()
coor = pdb_numpy.Coor("3eam.pdb")
pdb_numpy_read_time = time.time() - start_time 
print(f"Time taken to get coordinates:   {pdb_numpy_read_time:.4f} seconds")
start_time = time.time()
coor.write("tmp_python.pdb", overwrite=True)
pdb_numpy_write_time = time.time() - start_time
print(f"Time taken to write coordinates: {pdb_numpy_write_time:.4f} seconds")

start_time = time.time()
model = pdb_cpp.Model()
model.read("3eam.pdb")
cpp_read_time = time.time() - start_time
print(f"Time taken to get coordinates:   {cpp_read_time:.4f} seconds  {pdb_numpy_read_time/cpp_read_time:.2f}")
start_time = time.time()
model.write("tmp.pdb")
cpp_write_time = time.time() - start_time
print(f"Time taken to write coordinates: {cpp_write_time:.4f} seconds  {pdb_numpy_write_time/cpp_write_time:.2f}")
