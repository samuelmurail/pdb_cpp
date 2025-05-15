import pdb_numpy
import time


start_time = time.time()
coor = pdb_numpy.Coor("3eam.pdb")
end_time = time.time()
print(f"Time taken to get coordinates: {end_time - start_time} seconds")