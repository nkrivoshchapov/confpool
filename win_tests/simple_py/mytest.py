print("Hi")
import os
os.add_dll_directory("C:\\msys64\\mingw64\\bin")
os.add_dll_directory("C:\\msys64\\usr\\bin")
os.add_dll_directory(os.getcwd())
import confpool as m # pymolxyz

print("He gave me " + repr(m.give_five()))
print("Hi")
