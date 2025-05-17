import pdb_cpp
import time


start_time = time.time()
model = pdb_cpp.Model()
model.read("3eam.pdb")
end_time = time.time()
print(f"Time taken to get coordinates:   {end_time - start_time:.4f} seconds")
start_time = time.time()
