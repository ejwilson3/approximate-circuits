import sys
import matplotlib.pyplot as plt
from scipy.spatial import distance
import ast

# Some numpy and qiskit combinations send out a ridiculous number of deprecation warnings
# import warnings
# warnings.filterwarnings('ignore')

if len(sys.argv) < 2:
    print("Missing directory.")
    exit()

# Read the file to obtain the results.
syn_dir = sys.argv[1]
datafile = open(syn_dir + "/" + syn_dir + ".counts.txt", "r")
lines = datafile.read()
datafile.close()
lines = lines.split("\n")
# This holds the data for the approximate circuits; keys are circuit depth.
result_dict = ast.literal_eval(lines[0])
# This is the CNOT depth for the synthesis result circuit
fullcount = ast.literal_eval(lines[1])
# This is the data for the syntehsis result circuit
result_full = ast.literal_eval(lines[2])
# This is the CNOT depth for the reference circuit
refcount = ast.literal_eval(lines[3])
# This is the data for the reference circuit
result_ref = ast.literal_eval(lines[4])
# This is the data for the reference circuit run on a simulator without noise.
result_perf = ast.literal_eval(lines[5])
nq = 0
for key in result_perf:
    nq = len(key)
    break

# This begins the comparison stage. Take the counts from the noiseless simulation to create distributions to compare.

# Store the bins that are in the perfect results.
matches = []
# Store the results of the perfect circuit in a list usable by the Jensen Shannon calculation
target = []
# Store the Jensen Shannon distance of each circuti; key is the CNOT depth.
data_dict = {}
# This creates a list of possible results based on the number of qubits.
bins = ['0', '1']
for i in range(1, nq):
    bins2 = []
    for j, string in enumerate(bins):
        bins[j] = '0' + string
        bins2.append('1' + string)
    bins += bins2

# Build the matches and target lists
for key in result_perf:
    matches.append(key)
for targ in bins:
    if targ in matches:
        target.append(result_perf[targ])
    else:
        target.append(0.0)

# Actually compare the counts with the target
for key in result_dict:
    data_dict[key] = []
    for idx, step in enumerate(result_dict[key]):
        accumulated = []
        for match in bins:
            try:
                accumulated.append(step[match])
            except:
                accumulated.append(0.0)
        val = distance.jensenshannon(accumulated, target)
        data_dict[key].append(val)

accumulated = []
for match in bins:
    try:
        accumulated.append(result_full[match])
    except:
        accumulated.append(0.0)
data_full = distance.jensenshannon(accumulated, target)
    
accumulated = []   
for match in bins:
    try:
        accumulated.append(result_ref[match])
    except:
        accumulated.append(0.0)
data_ref = distance.jensenshannon(accumulated, target)


# Create the plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
xs = []
ys = []
for j, key in enumerate(data_dict):
    for k, val in enumerate(data_dict[key]):
        xs.append(key)
        ys.append(val)
if(len(xs) > 0):
    ax1.scatter(xs, ys, c='blue', s=10)

ax1.plot(range(max(xs) + 1), [data_ref]*(max(xs) + 1), lw=.5, linestyle=':')
ax1.scatter(fullcount, data_full, c='red', s=10)
ax1.scatter(refcount, data_ref, c='yellow', s=10)

ax1.set_ylabel("JS distance from expected value")
ax1.set_xlabel("CNOT count")
plt.savefig(syn_dir + "/" + syn_dir + ".png")
