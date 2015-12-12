import os
import subprocess
import math

names = []
nodeslist = []
pernode = 16

# Name of this run
run_name = "test2"


def enqueue_sim(name, sizeX, sizeY, sizeZ, nProcX, nProcY, nProcZ):
    # Read configuration file template
    inConf = open("conf_cavity_3D_unif_dns.xml", 'r')
    conf_file = inConf.read()
    inConf.close()

    # Substitute parameters
    conf_file = conf_file.format(run_name=run_name, name=name, sX=sizeX, sY=sizeY, sZ=sizeZ, npX=nProcX, npY=nProcY, npZ=nProcZ)

    # Write new configuration file
    outConf = open("tmp/{0}/{1}".format(run_name, name), 'w')
    outConf.write(conf_file)
    outConf.close()

    # Enqueue simulation
    names.append(name)
    nodeslist.append(nProcX * nProcY * nProcZ)

if __name__ == '__main__':

    # Create folders
    if not os.path.exists("tmp/{}".format(run_name)):
        os.mkdir("tmp/{}".format(run_name))
    if not os.path.exists("../../../output/{}".format(run_name)):
        os.mkdir("../../../output/{}".format(run_name))


    # Enqueue simulations to run:
    enqueue_sim("1", 20, 20, 20, 2, 2, 2)
    enqueue_sim("2", 20, 20, 20, 4, 4, 2)



    # Read batch file template
    inBatch = open("batch.sh", 'r')
    batch = inBatch.read()
    inBatch.close()

    # Substitute parameters
    batch = batch.format(nodes = int(math.ceil(max(nodeslist) / pernode)))
    # List all enqueued simulations
    for i in range(len(names)):
        batch += "mpirun -np {0} -ppn {1} ./ns conf/parallel/cavity/tmp/{2}/{3} 2>>output/{2}/times.txt\n".format(nodeslist[i], pernode, run_name, names[i])

    # Write new batch file
    outBatch = open("tmp/{}/batch.sh".format(run_name), 'w')
    outBatch.write(batch)
    outBatch.close()

    # Run all enqueued simulations
    subprocess.call("sbatch tmp/{}/batch.sh".format(run_name), shell=True)
