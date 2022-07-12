import os, sys, shutil, subprocess, glob

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

try:
    shutil.rmtree("release/")
except:
    pass
os.mkdir("release")
os.mkdir("release/myconfscript")
shutil.copy2("build/confpool.cpython-38-x86_64-linux-gnu.so", "release/myconfscript")
for file in glob.glob("release_files/myconfscript/*"):
    if os.path.isfile(file):
        shutil.copy2(file, "release/myconfscript/")
for file in glob.glob("release_files/*"):
    if os.path.isfile(file):
        shutil.copy2(file, "release/")
