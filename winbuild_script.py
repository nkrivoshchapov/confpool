import os, sys, shutil, subprocess, glob, ntpath

# raise Exception('E')
ver = int(sys.argv[1])
setup_lines = open("win_release_files/setup.py", 'r').readlines()
for i in range(len(setup_lines)):
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
subprocess.run('ninja', shell=True)
os.chdir(mainwd)

print("Creating release directory")
try:
    shutil.rmtree("release/")
except:
    pass
os.mkdir("release")
os.mkdir("release/pyxyz")
sofiles = glob.glob("build/confpool*.pyd")
assert len(sofiles) == 1
shutil.copy2(sofiles[0], "release/pyxyz/confpool.pyd")

shutil.copy2("release/pyxyz/confpool.pyd", "release/pyxyz/confpool.so")
mainwd = os.getcwd()
sopath = os.path.join(mainwd, "release/pyxyz/confpool.so")
sseq = ['ldd', sopath]
print(" ".join(sseq))
out = subprocess.run(sseq, stdout=subprocess.PIPE)
outlines = str(out.stdout).split("\\n")
deps = []
for line in outlines:
    # print(repr(line))
    if len(line) > 10:
        deps.append(line.split('=>')[1].split('(')[0].strip())

for dep in deps:
    print(dep)
    if not dep.startswith('/c/Windows/'):
        # assert os.path.isfile(dep)
        sseq = ['cp',dep, "release/pyxyz/"]
        subprocess.run(sseq, stdout=subprocess.PIPE)
os.remove("release/pyxyz/confpool.so")
for dllfile in glob.glob("save_dlls/*dll"):
    shutil.copy2(dllfile, "release/pyxyz/")

for file in glob.glob("win_release_files/pyxyz/*"):
    if os.path.isfile(file):
        shutil.copy2(file, "release/pyxyz/")
    elif os.path.isdir(file):
        shutil.copytree(file, "release/pyxyz/{}".format(ntpath.basename(file)))
for file in glob.glob("win_release_files/*"):
    if os.path.isfile(file):
        shutil.copy2(file, "release/")
print("Done creating release directory")

print("Updating setup.py")
with open("release/setup.py", 'w') as f:
    f.write(setup_text)
print("Done updating setup.py")
