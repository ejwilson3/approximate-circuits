import warnings
import qiskit
from qiskit import QuantumRegister
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from qiskit.providers.aer import noise
import re
from scipy.spatial import distance

# Currently, my numpy deprecation warnings with qiskit will eat up all the output
warnings.filterwarnings('ignore')

if len(sys.argv) < 4:
    print("Missing one of: reference QASM, synthesis directory, or number of synthesis files to read.")
    exit()

reference_name = sys.argv[1]
reference = qiskit.QuantumCircuit.from_qasm_file(sys.argv[1])
syn_dir = sys.argv[2]
n_files = int(sys.argv[3])

simulated = True
shots = 1000
# optimization level
ol = 1
qiskit.IBMQ.load_account()
provider = qiskit.IBMQ.get_provider(hub='ibm-q', group='open', project='main')
#gather fidelity statistics on this device if you want to create a noise model for the simulator
device = provider.get_backend('ibmq_5_yorktown')
properties = device.properties()
coupling_map = device.configuration().coupling_map
noise_model = noise.NoiseModel.from_backend(properties)
basis_gates = noise_model.basis_gates
simulator = qiskit.Aer.get_backend('qasm_simulator')


#For singular QASMs
base = syn_dir + "/" + syn_dir
tail = ".qasm"
g_circuits = []
# g_circuits_1 = []
# g_circuits_0 = []
g_depths = []
g_count = []
g_HS = []
maxdepth = 0
# n_files = 25
nq = reference.num_qubits

skipped = 0
for i in range(0, n_files+1):
    filename = base + str(i) + tail
    circ = qiskit.QuantumCircuit.from_qasm_file(filename)
    if(circ.num_qubits != nq):
        skipped += 1
        continue
    # matdat = qiskit.quantum_info.Operator(circ)
    # matrix = matdat.data
    # g_HS.append(1 - np.abs(np.sum(np.multiply(Orig,np.conj(matrix)))) / Orig.shape[0])
    # g_circuits_0.append(qiskit.QuantumCircuit(nq))
    # circ1 = qiskit.QuantumCircuit(nq)
    circ2 = qiskit.QuantumCircuit(nq)
#     for j in range(nq-1):
      #   if(j != 0):
        #     circ2.x(j)
 #        circ2.h(j)
    g_circuits.append(circ2)
    # g_circuits_1.append(circ1)
    g_circuits[i-skipped].extend(circ)
    # g_circuits_0[i].extend(circ)
    # g_circuits_1[i].extend(circ)
    cr = qiskit.ClassicalRegister(nq)
    g_circuits[i-skipped].add_register(cr)
    # g_circuits_0[i].add_register(cr)
    # g_circuits_1[i].add_register(cr)
    for j in range(nq):
        g_circuits[i-skipped].measure(j, cr[j])
        # g_circuits_0[i].measure(j, cr[j])
        # g_circuits_1[i].measure(j, cr[j])

g_full_circuit = g_circuits.pop()

# g_nat_circuit = qiskit.QuantumCircuit(nq)
# g_nat_circuit_1 = qiskit.QuantumCircuit(nq)
# g_nat_circuit_0 = qiskit.QuantumCircuit(nq)
# for j in range(nq-1):
    # if(j != 0):
        # g_nat_circuit.x(j)
    # g_nat_circuit_1.x(j)
g_nat_circuit = reference
# g_nat_circuit_0.extend(circ)
# g_nat_circuit_1.extend(circ)
cr = qiskit.ClassicalRegister(nq)
g_nat_circuit.add_register(cr)
# g_nat_circuit_0.add_register(cr)
# g_nat_circuit_1.add_register(cr)
for j in range(nq):
    g_nat_circuit.measure(j, cr[j])
    # g_nat_circuit_0.measure(j, cr[j])
    # g_nat_circuit_1.measure(j, cr[j])

    
for idx, circ in enumerate(g_circuits):
    tcirc = qiskit.compiler.transpile(circ, backend=device, basis_gates=basis_gates, optimization_level=ol)
    qobj = qiskit.compiler.assemble(tcirc, shots=shots)
    g_circuits[idx] = qobj

g_count = []
for circ in g_circuits:
    n_cx = 0
    for gate in circ.to_dict()['experiments'][0]['instructions']:
        if(len(gate['qubits']) > 1):
            n_cx += 1
    g_count.append(n_cx)

tcirc = qiskit.compiler.transpile(g_full_circuit, backend=device, optimization_level=ol, basis_gates=basis_gates)
qobj1 = qiskit.compiler.assemble(tcirc, shots=shots)

tcirc = qiskit.compiler.transpile(g_nat_circuit, backend=device, optimization_level=ol, basis_gates=basis_gates)
qobj2 = qiskit.compiler.assemble(tcirc, shots=shots)

g_completed = {}
# g_completed_0 = {}
# g_completed_1 = {}
# print(circuits)
# for key in g_raw_counts:
    # g_completed[key] = []
for key in range(min(g_count), max(g_count) + 1):
    g_completed[key] = []
  #   g_completed_0[key] = []
    # g_completed_1[key] = []
    
for i, depth in enumerate(g_count):
    if simulated:
        g_completed[depth].append(simulator.run(g_circuits[i], noise_model=noise_model))
        # g_completed_0[depth].append(qiskit.execute(g_circuits_0[i], simulator, noise_model=noise_model, coupling_map=coupling_map,basis_gates=basis_gates))
      #   g_completed_0[depth].append(simulator.run(g_circuits_0[i], noise_model=noise_model))
      #   g_completed_1[depth].append(simulator.run(g_circuits_1[i], noise_model=noise_model))
    else:
        g_completed[depth].append(device.run(g_circuits[i]))

if simulated:
    g_full_completed = simulator.run(qobj1, noise_model=noise_model)
    g_nat_completed = simulator.run(qobj2, noise_model=noise_model)
else:
    g_full_completed = device.run(qobj1)
    g_nat_completed = device.run(qobj2)

fullcount = 0
natcount = 0
for gate in qobj1.to_dict()['experiments'][0]['instructions']:
    if(len(gate['qubits']) > 1):
        fullcount += 1

for gate in qobj2.to_dict()['experiments'][0]['instructions']:
    if(len(gate['qubits']) > 1):
        natcount += 1

perfect = simulator.run(qobj2)

g_anotherone = {}
# g_anotherone_0 = {}
# g_anotherone_1 = {}
for key in g_completed:
    g_anotherone[key] = []
  #   g_anotherone_0[key] = []
   #  g_anotherone_1[key] = []
for key in g_completed:
    for res in g_completed[key]:
        if res != []:
            g_anotherone[key].append(res.result().get_counts())
   #  for res in g_completed_0[key]:
     #    if res != []:
       #      g_anotherone_0[key].append(res.result().get_counts())
    # for res in g_completed_1[key]:
      #   if res != []:
        #     g_anotherone_1[key].append(res.result().get_counts())


g_another_full = g_full_completed.result().get_counts()
g_another_nat = g_nat_completed.result().get_counts()
perfect_res = perfect.result().get_counts()

bins = ['0', '1']
for i in range(1, nq):
    bins2 = []
    for j, string in enumerate(bins):
        bins[j] = '0' + string
        bins2.append('1' + string)
    bins += bins2
matches = []
target = []
for key in perfect_res:
    matches.append(key)
for targ in bins:
    if targ in matches:
        target.append(perfect_res[targ])
    else:
        target.append(0.0)

g_data = {}
for key in g_anotherone:
    # print(key)
    g_data[key] = []
    i = 0
    j = 0
    for idx, step in enumerate(g_anotherone[key]):
        # temp = []
        # print(step)
        # for k, circ in enumerate(step):
            # print(circ)
            # if(type(circ) == type([])):
                # print("HA")
                # for piece in circ:
                    # temp.append(piece['010']/1024)
            # else:
        accumulated = []
        for match in bins:
            try:
            # temp.append(step['1011']/1024)
                accumulated.append(step[match])
            
            # temp.append(step['000']/1000)
            except:
                accumulated.append(0.0)
        # temp.append(accumulated/1000.0)
        # print(temp)
        val = distance.jensenshannon(accumulated, target)
        g_data[key].append(val)

accumulated = []
for match in bins:
    try:
    # g_data_full = g_another_full['1011']/1024
        # accumulated = g_another_full['0010'] + g_another_full['0000'] + g_another_full['0001'] + \
          #                                      g_another_full['1111'] + g_another_full['0100'] + g_another_full['0101'] + g_another_full['0110'] + g_another_full['0011']
        accumulated.append(g_another_full[match])
    except:
        accumulated.append(0.0)
    
g_data_full = distance.jensenshannon(accumulated, target)
    # g_data_full = g_another_full['000']/1000
    
accumulated = []   
for match in bins:
    try:
    # g_data_nat = g_another_nat['000 000']/1000
    # accumulated = g_another_nat['0010'] + g_another_nat['0000'] + g_another_nat['1000'] + \
      #                                         g_another_nat['1111'] + g_another_nat['0100'] + g_another_nat['0101'] + g_another_nat['0110'] + g_another_nat['0011']
        accumulated.append(g_another_nat[match])
    except:
        accumulated.append(0.0)
g_data_nat = distance.jensenshannon(accumulated, target)

# duplicates = []
# for job in jobs:
    # counts = job.result().get_counts()
    # if type(counts) == type({}):
        # duplicates.append(1)
    # else:
        # duplicates.append(len(counts))

i = 0
fig = plt.figure()
ax1 = fig.add_subplot(111)
xs = []
ys = []
for j, key in enumerate(g_data):
# for j, key in enumerate(data):
    for k, val in enumerate(g_data[key]):
        xs.append(key)
        ys.append(val)
        # xs.append((j+1))#*delta_t)
        # ys.append(avg_mag_qc[0][i])
        # i += 1
if(len(xs) > 0):
    ax1.scatter(xs, ys, c='blue', s=10)

ax1.plot(range(max(xs) + 1), [g_data_nat]*(max(xs) + 1), lw=.5, linestyle=':')
ax1.scatter(fullcount, g_data_full, c='red', s=10)
ax1.scatter(natcount, g_data_nat, c='yellow', s=10)

ax1.set_ylabel("JS distance from expected value")
ax1.set_xlabel("CNOT count")
plt.savefig(syn_dir + "/" + syn_dir + ".png")

if not simulated:
    j_dict = {}
    for key in g_completed:
        j_dict[key] = []
        for item in g_completed[key]:
            j_dict[key].append(item.job_id())
    newfile = open(syn_dir + "/" + syn_dir + ".backupJobIDs.txt", "w")
    newfile.write(str(j_dict) + "\n" + str(g_full_completed.job_id()) + "\n" + str(g_nat_completed.job_id()) + "\n")
    newfile.close()
