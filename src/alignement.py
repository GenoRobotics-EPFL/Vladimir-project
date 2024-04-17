

from src.dataRetriever import *

reads = getReads(11)

pathInput = os.path.join("fastqpass", "iteration11.fastq")
pathResult = os.path.join("tmp_result.paf")

minimap2_command = f"minimap2 -x ava-ont {pathInput} {pathInput} > {pathResult}"
os.system(minimap2_command)

# * find the read with the highest number of mappings

# readName -> [lines]
occurences = {}

with open(pathResult, "r") as result:


    while True:
        line = result.readline()

        if line == "":
            break

        line = line[:-1] # remove newline

        objects = line.split('\t')

        nameRead = objects[0]

        if nameRead in occurences:
            occurences[nameRead].append(objects)
        else:
            occurences[nameRead] = [objects]

readHighestMapping = ""
highestMapping = 0

for k, v in occurences.items():
    if len(v) > highestMapping:
        highestMapping = len(v)
        readHighestMapping = k

# visualise the read and the first mapping





