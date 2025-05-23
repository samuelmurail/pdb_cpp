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
select_times = []
read_cpp_times = []
write_cpp_times = []
select_cpp_times = []

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
    start_time = time.time()
    selection = "resid > 250 and chain A B and resname ALA GLY and within 20 of chain C"
    coor.select_atoms(selection)
    pdb_numpy_select_time = time.time() - start_time
    select_times.append(pdb_numpy_select_time)

    # cpp
    start_time = time.time()
    coor = core.Coor(file_name)
    pdb_numpy_read_time = time.time() - start_time 
    read_cpp_times.append(pdb_numpy_read_time)
    start_time = time.time()
    coor.write(f"tmp_cpp_{pdb_id}.pdb")
    pdb_numpy_write_time = time.time() - start_time
    write_cpp_times.append(pdb_numpy_write_time)

    start_time = time.time()
    coor.select_atoms(selection)
    pdb_numpy_select_time = time.time() - start_time
    select_cpp_times.append(pdb_numpy_select_time)


avg_read, std_read = avg_std(read_times)
print(f"-pdb_numpy Time taken to get coordinates:   {avg_read:.4f} +- {std_read:.4f} seconds")
avg_write, std_write = avg_std(write_times)
print(f"-pdb_numpy Time taken to write coordinates: {avg_write:.4f} +- {std_write:.4f} seconds")
avg_select, std_select = avg_std(select_times)
print(f"-pdb_numpy Time taken to select atoms:      {avg_select:.4f} +- {std_select:.4f} seconds")
avg_cpp_read, std_cpp_read = avg_std(read_cpp_times)
print(f"-pdb_cpp   Time taken to get coordinates:   {avg_cpp_read:.4f} +- {std_cpp_read:.4f} seconds, speed-up:  {avg_read/avg_cpp_read:.2f}")
avg_cpp_write, std_cpp_write = avg_std(write_cpp_times)
print(f"-pdb_cpp   Time taken to write coordinates: {avg_cpp_write:.4f} +- {std_cpp_write:.4f} seconds, speed-up:  {avg_write/avg_cpp_write:.2f} ")
avg_cpp_select, std_select = avg_std(select_cpp_times)
print(f"-pdb_cpp   Time taken to select atoms:      {avg_cpp_select:.4f} +- {std_select:.4f} seconds, speed-up:  {avg_select/avg_cpp_select:.2f}")

read_times = []
write_times = []

file_name = "2rri.pdb"
pdb_id = file_name.split(".")[0]
read_times = []
write_times = []
select_times = []
read_cpp_times = []
write_cpp_times = []
select_cpp_times = []

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
    start_time = time.time()
    selection = "resid <= 50 and chain A and resname ALA GLY CYS"
    coor.select_atoms(selection)
    pdb_numpy_select_time = time.time() - start_time
    select_times.append(pdb_numpy_select_time)

    # cpp
    start_time = time.time()
    coor = core.Coor(file_name)
    pdb_numpy_read_time = time.time() - start_time 
    read_cpp_times.append(pdb_numpy_read_time)
    start_time = time.time()
    coor.write(f"tmp_cpp_{pdb_id}.pdb")
    pdb_numpy_write_time = time.time() - start_time
    write_cpp_times.append(pdb_numpy_write_time)

    start_time = time.time()
    coor.select_atoms(selection)
    pdb_numpy_select_time = time.time() - start_time
    select_cpp_times.append(pdb_numpy_select_time)


avg_read, std_read = avg_std(read_times)
print(f"-pdb_numpy Time taken to get coordinates:   {avg_read:.4f} +- {std_read:.4f} seconds")
avg_write, std_write = avg_std(write_times)
print(f"-pdb_numpy Time taken to write coordinates: {avg_write:.4f} +- {std_write:.4f} seconds")
avg_select, std_select = avg_std(select_times)
print(f"-pdb_numpy Time taken to select atoms:      {avg_select:.4f} +- {std_select:.4f} seconds")
avg_cpp_read, std_cpp_read = avg_std(read_cpp_times)
print(f"-pdb_cpp   Time taken to get coordinates:   {avg_cpp_read:.4f} +- {std_cpp_read:.4f} seconds, speed-up:  {avg_read/avg_cpp_read:.2f}")
avg_cpp_write, std_cpp_write = avg_std(write_cpp_times)
print(f"-pdb_cpp   Time taken to write coordinates: {avg_cpp_write:.4f} +- {std_cpp_write:.4f} seconds, speed-up:  {avg_write/avg_cpp_write:.2f} ")
avg_cpp_select, std_select = avg_std(select_cpp_times)
print(f"-pdb_cpp   Time taken to select atoms:      {avg_cpp_select:.4f} +- {std_select:.4f} seconds, speed-up:  {avg_select/avg_cpp_select:.2f}")

print("3eam:")
file_name = "3eam.pdb"
coor = pdb_numpy.Coor(file_name)
# print(coor.transformation)
# print(coor.symmetry)

sel = coor.select_atoms("name CA CB CG")
print(F"Number of CA CB CG : {sel.len}")
sel = coor.select_atoms("resname ALA GLY CYS")
print(F"Number of ALA GLY CYS: {sel.len}")

sel = coor.select_atoms("x >= 30")
print(F"Number of x>=30 : {sel.len}")
sel = coor.select_atoms("z < 20")
print(F"Number of z <20 : {sel.len}")
sel = coor.select_atoms("resid < 20.0")
print(F"Number of resid <20 : {sel.len}")

sel = coor.select_atoms("resid 20 21 25")
print(F"Number of resid 20 21 25 : {sel.len}")


sel = coor.select_atoms("chain A B")
print(F"Number of chain A B : {sel.len}")

sel = coor.select_atoms("name C*")
print(F"Number of name C* : {sel.len}")

sel = coor.select_atoms("resid > 250 and chain A B and resname ALA GLY")
print(F"Number of :{selection} : {sel.len}")

selection = "resid > 250 and chain A B and resname ALA GLY and within 20.0 of chain C";
sel = coor.select_atoms(selection)
print(F"Number of :{selection} : {sel.len}")

selection = "resid >= 250 and not chain C D E and resname ALA GLY and within 20.0 of chain C"
sel = coor.select_atoms(selection)
print(F"Number of :{selection} : {sel.len}")


selection = "within 10.0 of chain C";
sel = coor.select_atoms(selection)
print(F" within 10.0 of chain C : {sel.len}")

print("3rri:")
file_name = "2rri.pdb"
coor = pdb_numpy.Coor(file_name)
# print(coor.transformation)
# print(coor.symmetry)


#print(coor.models[0].chain)

coor = core.Coor(file_name)
#print(coor.get_Models(0).get_chain())

selection = "resid >= 250 and not chain C D E and resname ALA GLY and within 20.0 of chain C"
from pdb_numpy import select

print(select.parse_selection(selection))

