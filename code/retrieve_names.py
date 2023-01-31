import os

# getting the cwd brings up the Snakefile path instead of the script path, then the rest of the path
# is taken from "params" in the snakemake file, everything is joined together in strings when needed

script_path = os.getcwd()
in_dir = snakemake.params.trimmed_files_loc

os.chdir(script_path + "/" + in_dir)
files = os.listdir()

seq = list()

for i in files:
    if ".gz" in i:
        seq.append(i)

R1 = list()

for i in seq:
    if "R1.f" in i:
        R1.append(i)

names = list()

for i in R1:
    name = i.split(".")[0]
    names.append(name)

file = open('samples.txt', 'w')
new_content = '\n'.join(names)
file.write(new_content)
locals()
file.close()

exit()
