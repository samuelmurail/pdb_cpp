import time
import pdb_numpy
import sys

sys.path.insert(0, "./src")

from pdb_cpp import core


file_name = "3eam.pdb"

start_time = time.time()
coor = pdb_numpy.Coor(file_name)
pdb_numpy_read_time = time.time() - start_time 
print(f"Time taken to get coordinates:   {pdb_numpy_read_time:.4f} seconds")
start_time = time.time()
coor.write("tmp_python.pdb", overwrite=True)
pdb_numpy_write_time = time.time() - start_time
print(f"Time taken to write coordinates: {pdb_numpy_write_time:.4f} seconds")

start_time = time.time()
coor = core.Coor()
coor.read(file_name)
cpp_read_time = time.time() - start_time
print(f"Time taken to get coordinates:   {cpp_read_time:.4f} seconds  {pdb_numpy_read_time/cpp_read_time:.2f}")
start_time = time.time()
coor.write("tmp.pdb")
cpp_write_time = time.time() - start_time
print(f"Time taken to write coordinates: {cpp_write_time:.4f} seconds  {pdb_numpy_write_time/cpp_write_time:.2f}")



file_name = "2rri.pdb"
start_time = time.time()
coor = pdb_numpy.Coor(file_name)
pdb_numpy_read_time = time.time() - start_time 
print(f"Time taken to get coordinates:   {pdb_numpy_read_time:.4f} seconds")
start_time = time.time()
coor.write("tmp_python.pdb", overwrite=True)
pdb_numpy_write_time = time.time() - start_time
print(f"Time taken to write coordinates: {pdb_numpy_write_time:.4f} seconds")

start_time = time.time()
coor = core.Coor()
coor.read(file_name)
cpp_read_time = time.time() - start_time
print(f"Time taken to get coordinates:   {cpp_read_time:.4f} seconds  {pdb_numpy_read_time/cpp_read_time:.2f}")
start_time = time.time()
coor.write("tmp.pdb")
cpp_write_time = time.time() - start_time
print(f"Time taken to write coordinates: {cpp_write_time:.4f} seconds  {pdb_numpy_write_time/cpp_write_time:.2f}")
