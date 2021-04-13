import qiskit
# import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from qiskit.providers.aer import noise
# import re
from scipy.spatial import distance

# Some numpy and qiskit combinations send out a ridiculous number of deprecation warnings
# import warnings
# warnings.filterwarnings('ignore')

if len(sys.argv) < 4:
    print("Missing one of: reference QASM, synthesis directory, or number of synthesis files to read.")
    exit()

################################################################################
# Most editing will need to be done here
simulated = True
shots = 1000
# optimization level
ol = 1
qiskit.IBMQ.load_account()
provider = qiskit.IBMQ.get_provider(hub='ibm-q', group='open', project='main')
device = provider.get_backend('ibmq_5_yorktown')
################################################################################

reference = qiskit.QuantumCircuit.from_qasm_file(sys.argv[1])
syn_dir = sys.argv[2]
n_files = int(sys.argv[3])
nq = reference.num_qubits

################################################################################
# If the reference circuit is measured, comment this section.
cr = qiskit.ClassicalRegister(nq)
reference.add_register(cr)
for j in range(nq):
    reference.measure(j, cr[j])
################################################################################

properties = device.properties()
coupling_map = device.configuration().coupling_map
noise_model = noise.NoiseModel.from_backend(properties)
basis_gates = noise_model.basis_gates
simulator = qiskit.Aer.get_backend('qasm_simulator')

# Generate the file names
base = syn_dir + "/" + syn_dir
tail = ".qasm"

# Store the circuits themselves
g_circuits = []
# Store the CNOT depth of each circuit
g_count = []
fullcount = 0
natcount = 0
# Store the completed circuits; the keys are the CNOT depth.
g_completed = {}
# Store the results of the circuits; the keys are the CNOT depth.
result_dict = {}

skipped = 0
for i in range(0, n_files+1):
    filename = base + str(i) + tail
    circ = qiskit.QuantumCircuit.from_qasm_file(filename)
    if(circ.num_qubits != nq):
        # This circuit has too few (or too many) qubits and won't be comparable by most metrics. QFast will sometimes create smaller circuits.
        skipped += 1
        continue
    circ2 = qiskit.QuantumCircuit(nq)
    g_circuits.append(circ2)
    g_circuits[i-skipped].extend(circ)
    # This circuits from QSearch and QFast don't come with measurement
    cr = qiskit.ClassicalRegister(nq)
    g_circuits[i-skipped].add_register(cr)
    for j in range(nq):
        g_circuits[i-skipped].measure(j, cr[j])

# This is the circuit normally output by the synthesizer.
g_full_circuit = g_circuits.pop()

# Transpile the approximate circuits. It may fail here, if the transpilation
# times out. You can try again, or lower the optimization level if this happens
# too often.    
for idx, circ in enumerate(g_circuits):
    tcirc = qiskit.compiler.transpile(circ, backend=device, basis_gates=basis_gates, optimization_level=ol)
    qobj = qiskit.compiler.assemble(tcirc, shots=shots)
    g_circuits[idx] = qobj

tcirc = qiskit.compiler.transpile(g_full_circuit, backend=device, optimization_level=ol, basis_gates=basis_gates)
qobj1 = qiskit.compiler.assemble(tcirc, shots=shots)

tcirc = qiskit.compiler.transpile(reference, backend=device, optimization_level=ol, basis_gates=basis_gates)
qobj2 = qiskit.compiler.assemble(tcirc, shots=shots)

# Actually counting the number of CNOTs; transpilation can add some.
for circ in g_circuits:
    n_cx = 0
    for gate in circ.to_dict()['experiments'][0]['instructions']:
        if(len(gate['qubits']) > 1):
            n_cx += 1
    g_count.append(n_cx)

for gate in qobj1.to_dict()['experiments'][0]['instructions']:
    if(len(gate['qubits']) > 1):
        fullcount += 1

for gate in qobj2.to_dict()['experiments'][0]['instructions']:
    if(len(gate['qubits']) > 1):
        natcount += 1


for key in range(min(g_count), max(g_count) + 1):
    g_completed[key] = []
    
# Here actually run the circuits
if simulated:
    for i, depth in enumerate(g_count):
        g_completed[depth].append(simulator.run(g_circuits[i], noise_model=noise_model))
    g_full_completed = simulator.run(qobj1, noise_model=noise_model)
    g_nat_completed = simulator.run(qobj2, noise_model=noise_model)
else:
    for i, depth in enumerate(g_count):
        g_completed[depth].append(device.run(g_circuits[i]))
    g_full_completed = device.run(qobj1)
    g_nat_completed = device.run(qobj2)

# This will be run on the simulator to compare. This is how the JS score is generated; it is compared to the results of this circuit run without noise.
perfect = simulator.run(qobj2)

# Put the counts in their own dictionaries
for key in g_completed:
    result_dict[key] = []
for key in g_completed:
    for res in g_completed[key]:
        if res != []:
            result_dict[key].append(res.result().get_counts())
result_full = g_full_completed.result().get_counts()
result_nat = g_nat_completed.result().get_counts()
result_perf = perfect.result().get_counts()


# This begins the comparison stage. Take the counts from the noiseless simulation to create distributions to compare.
# Store the bins that are in the perfect results.
matches = []
# Store the results of the perfect circuit in a list usable by the Jensen Shannon calculation
target = []
# Store the Jensen Shannon distance of each circuti; key is the CNOT depth.
g_data = {}
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
    g_data[key] = []
    for idx, step in enumerate(result_dict[key]):
        accumulated = []
        for match in bins:
            try:
                accumulated.append(step[match])
            except:
                accumulated.append(0.0)
        val = distance.jensenshannon(accumulated, target)
        g_data[key].append(val)

accumulated = []
for match in bins:
    try:
        accumulated.append(result_full[match])
    except:
        accumulated.append(0.0)
g_data_full = distance.jensenshannon(accumulated, target)
    
accumulated = []   
for match in bins:
    try:
        accumulated.append(result_nat[match])
    except:
        accumulated.append(0.0)
g_data_nat = distance.jensenshannon(accumulated, target)


# Create the plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
xs = []
ys = []
for j, key in enumerate(g_data):
    for k, val in enumerate(g_data[key]):
        xs.append(key)
        ys.append(val)
if(len(xs) > 0):
    ax1.scatter(xs, ys, c='blue', s=10)

ax1.plot(range(max(xs) + 1), [g_data_nat]*(max(xs) + 1), lw=.5, linestyle=':')
ax1.scatter(fullcount, g_data_full, c='red', s=10)
ax1.scatter(natcount, g_data_nat, c='yellow', s=10)

ax1.set_ylabel("JS distance from expected value")
ax1.set_xlabel("CNOT count")
plt.savefig(syn_dir + "/" + syn_dir + ".png")


# Store the job IDs if run on a physical machine.
if not simulated:
    j_dict = {}
    for key in g_completed:
        j_dict[key] = []
        for item in g_completed[key]:
            j_dict[key].append(item.job_id())
    newfile = open(syn_dir + "/" + syn_dir + ".backupJobIDs.txt", "w")
    newfile.write(str(j_dict) + "\n" + str(g_full_completed.job_id()) + "\n" + str(g_nat_completed.job_id()) + "\n")
    newfile.close()
