import os, sys, shutil, subprocess, glob, ntpath

# name_ver = sys.argv[1] # like 'py37'
ver = int(sys.argv[1]) # '7' in this case
setup_lines = open("release_files/setup.py", 'r').readlines()
for i in range(len(setup_lines)):
    # if '{name_ver}' in setup_lines[i]:
    #     setup_lines[i] = setup_lines[i].format(name_ver=name_ver)
    if '{ver}' in setup_lines[i]:
        setup_lines[i] = setup_lines[i].format(ver=ver)
setup_text = ''.join(setup_lines)

try:
    shutil.rmtree("build/")
except:
    pass
os.mkdir("build")
mainwd = os.getcwd()
os.chdir("./build")
subprocess.run('cmake ..', shell=True)
subprocess.run('make', shell=True)
os.chdir(mainwd)

print("Creating release directory")
try:
    shutil.rmtree("release/")
except:
    pass
os.mkdir("release")
os.mkdir("release/pyxyz")
sofiles = glob.glob("build/confpool.*.so")
assert len(sofiles) == 1
shutil.copy2(sofiles[0], "release/pyxyz/confpool.so")
for file in glob.glob("release_files/pyxyz/*"):
    if os.path.isfile(file):
        shutil.copy2(file, "release/pyxyz/")
    elif os.path.isdir(file):
        shutil.copytree(file, "release/pyxyz/{}".format(ntpath.basename(file)))
for file in glob.glob("release_files/*"):
    if os.path.isfile(file):
        shutil.copy2(file, "release/")
print("Done creating release directory")

print("Updating setup.py")
with open("release/setup.py", 'w') as f:
    f.write(setup_text)
print("Done updating setup.py")
