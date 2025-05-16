import pdb_numpy
import time


start_time = time.time()
coor = pdb_numpy.Coor("3eam.pdb")
end_time = time.time()
print(f"Time taken to get coordinates:   {end_time - start_time:.4f} seconds")
start_time = time.time()
coor.write("tmp_python.pdb", overwrite=True)
end_time = time.time()
print(f"Time taken to write coordinates: {end_time - start_time:.4f} seconds")