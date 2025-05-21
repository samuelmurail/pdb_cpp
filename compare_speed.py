import time
import pdb_numpy
import sys

sys.path.insert(0, "./src")

from pdb_cpp import core

def avg_std(arr):
    """
    Calculate the average and standard deviation of a list of numbers.
    """
    avg = sum(arr) / len(arr)
    std = (sum((x - avg) ** 2 for x in arr) / len(arr)) ** 0.5
    return avg, std

N = 10

file_name = "3eam.pdb"
pdb_id = file_name.split(".")[0]
read_times = []
write_times = []
read_cpp_times = []
write_cpp_times = []

print(f"- Testing with {file_name} file")

for i in range(N):

    # pdb_numpy
    start_time = time.time()
    coor = pdb_numpy.Coor(file_name)
    pdb_numpy_read_time = time.time() - start_time 
    read_times.append(pdb_numpy_read_time)
    start_time = time.time()
    coor.write(f"tmp_python_{pdb_id}.pdb", overwrite=True)
    pdb_numpy_write_time = time.time() - start_time
    write_times.append(pdb_numpy_write_time)

    # cpp
    start_time = time.time()
    coor = core.Coor(file_name)
    pdb_numpy_read_time = time.time() - start_time 
    read_cpp_times.append(pdb_numpy_read_time)
    start_time = time.time()
    coor.write(f"tmp_cpp_{pdb_id}.pdb")
    pdb_numpy_write_time = time.time() - start_time
    write_cpp_times.append(pdb_numpy_write_time)


avg_read, std_read = avg_std(read_times)
print(f"-pdb_numpy Time taken to get coordinates:   {avg_read:.4f} +- {std_read:.4f} seconds")
avg_write, std_write = avg_std(write_times)
print(f"-pdb_numpy Time taken to write coordinates: {avg_write:.4f} +- {std_write:.4f} seconds")
avg_cpp_read, std_cpp_read = avg_std(read_cpp_times)
print(f"-pdb_cpp   Time taken to get coordinates:   {avg_cpp_read:.4f} +- {std_cpp_read:.4f} seconds, speed-up:  {avg_read/avg_cpp_read:.2f}")
avg_cpp_write, std_cpp_write = avg_std(write_cpp_times)
print(f"-pdb_cpp   Time taken to write coordinates: {avg_cpp_write:.4f} +- {std_cpp_write:.4f} seconds, speed-up:  {avg_write/avg_cpp_write:.2f} ")

read_times = []
write_times = []

file_name = "2rri.pdb"
pdb_id = file_name.split(".")[0]
read_times = []
write_times = []
read_cpp_times = []
write_cpp_times = []

print(f"- Testing with {file_name} file")

for i in range(N):

    # pdb_numpy
    start_time = time.time()
    coor = pdb_numpy.Coor(file_name)
    pdb_numpy_read_time = time.time() - start_time 
    read_times.append(pdb_numpy_read_time)
    start_time = time.time()
    coor.write(f"tmp_python_{pdb_id}.pdb", overwrite=True)
    pdb_numpy_write_time = time.time() - start_time
    write_times.append(pdb_numpy_write_time)

    # cpp
    start_time = time.time()
    coor = core.Coor(file_name)
    pdb_numpy_read_time = time.time() - start_time 
    read_cpp_times.append(pdb_numpy_read_time)
    start_time = time.time()
    coor.write(f"tmp_cpp_{pdb_id}.pdb")
    pdb_numpy_write_time = time.time() - start_time
    write_cpp_times.append(pdb_numpy_write_time)


avg_read, std_read = avg_std(read_times)
print(f"-pdb_numpy Time taken to get coordinates:   {avg_read:.4f} +- {std_read:.4f} seconds")
avg_write, std_write = avg_std(write_times)
print(f"-pdb_numpy Time taken to write coordinates: {avg_write:.4f} +- {std_write:.4f} seconds")
avg_cpp_read, std_cpp_read = avg_std(read_cpp_times)
print(f"-pdb_cpp   Time taken to get coordinates:   {avg_cpp_read:.4f} +- {std_cpp_read:.4f} seconds, speed-up:  {avg_read/avg_cpp_read:.2f}")
avg_cpp_write, std_cpp_write = avg_std(write_cpp_times)
print(f"-pdb_cpp   Time taken to write coordinates: {avg_cpp_write:.4f} +- {std_cpp_write:.4f} seconds, speed-up:  {avg_write/avg_cpp_write:.2f} ")

print("3eam:")
file_name = "3eam.pdb"
coor = pdb_numpy.Coor(file_name)
print(coor.transformation)
print(coor.symmetry)


print("3rri:")
file_name = "2rri.pdb"
coor = pdb_numpy.Coor(file_name)
print(coor.transformation)
print(coor.symmetry)

