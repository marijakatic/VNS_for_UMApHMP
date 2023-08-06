import subprocess

ns = [10, 15, 20, 25, 40, 50, 100, 200]
ps = [2, 3, 4, 5, 8]

for n in ns:
    for p in ps:
        subprocess.call(["./generate", str(n), str(p), "./generated/" + str(n) + "." + str(p)])