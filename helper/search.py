import os

for file in os.listdir("/media/data2/3d-dna-master"):
    script=False
    if os.path.isdir(f"/media/data2/3d-dna-master/{file}"):
        for subfile in os.listdir(f"/media/data2/3d-dna-master/{file}"):
            if subfile.endswith(".sh"):
                with open(f"/media/data2/3d-dna-master/{file}/{subfile}","r") as rfile:
                    line=rfile.read()
                    if "norms" in line:
                        print(subfile)
    elif file.endswith(".sh"):
        with open(f"/media/data2/3d-dna-master/{file}","r") as rfile:
            line=rfile.read()
            if "norms" in line:
                print(file)
