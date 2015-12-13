import os
import subprocess
import math

names = []
nodeslist = []
pernode = 16


def enqueue_sim(conf_template, name, solver, mesh, dim, lenX, lenY, lenZ, sizeX, sizeY, sizeZ, nProcX, nProcY, nProcZ):
    # Read configuration file template
    inConf = open(conf_template, 'r')
    conf_file = inConf.read()
    inConf.close()

    # Substitute parameters
    conf_file = conf_file.format(run_name=run_name, name=name, solver=solver, mesh=mesh, dim=dim, lX=lenX, lY=lenY, lZ=lenZ, sX=sizeX, sY=sizeY, sZ=sizeZ, npX=nProcX, npY=nProcY, npZ=nProcZ)

    # Write new configuration file
    outConf = open("tmp/{0}/{1}.xml".format(run_name, name), 'w')
    outConf.write(conf_file)
    outConf.close()

    # Enqueue simulation
    names.append(name)
    nodeslist.append(nProcX * nProcY * nProcZ)

def scaling():
    sX = 40
    sY = 40
    sZ = 40
    lX = 1.0
    lY = 1.0
    lZ = 1.0
    enqueue_sim(cav, "1", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "2_x", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 1, 1)
    enqueue_sim(cav, "2_y", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 2, 1)
    enqueue_sim(cav, "2_z", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 2)
    enqueue_sim(cav, "4", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 1)
    enqueue_sim(cav, "8", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)
    enqueue_sim(cav, "16_x", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)
    enqueue_sim(cav, "16_y", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 4, 2)
    enqueue_sim(cav, "16_z", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 4)
    enqueue_sim(cav, "32", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 4, 2)
    enqueue_sim(cav, "64", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 4, 4)

def validation():
    sX = 20
    sY = 20
    sZ = 20
    lX = 1.0
    lY = 1.0
    lZ = 1.0
    # enqueue_sim(cav, "cav_2D_seq", "dns", "uniform", 2, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    # enqueue_sim(cav, "cav_2D", "dns", "uniform", 2, lX, lY, lZ, sX, sY, sZ, 4, 2, 1)
    # enqueue_sim(cav, "cav_3D_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    # enqueue_sim(cav, "cav_3D_2", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 1, 1)
    # enqueue_sim(cav, "cav_3D_8", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)
    # enqueue_sim(cav, "cav_3D_16", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)
    lY = 5.0
    enqueue_sim(cav, "cav_y5.0_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "cav_y5.0", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)
    lY = 1.0
    lZ = 0.5
    enqueue_sim(cav, "cav_z0.5_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "cav_z0.5", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)

    sX = 10
    sY = 20
    sZ = 10
    lX = 5.0
    lY = 1.0
    lZ = 1.0
    enqueue_sim(bfs, "bfs_unif_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(bfs, "bfs_unif", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)
    enqueue_sim(bfs, "bfs_str_seq", "dns", "stretched", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(bfs, "bfs_str", "dns", "stretched", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)

if __name__ == '__main__':

    
    # Name of this run
    run_name = "validation_1"


    # Create folders
    if not os.path.exists("tmp/{}".format(run_name)):
        os.mkdir("tmp/{}".format(run_name))
    if not os.path.exists("../../output/{}".format(run_name)):
        os.mkdir("../../output/{}".format(run_name))


    # Enqueue simulations to run
    cav = "conf_cavity.xml"
    bfs = "conf_bfs.xml"
    # scaling()
    validation()


    # Read batch file template
    inBatch = open("batch.sh", 'r')
    batch = inBatch.read()
    inBatch.close()

    # Substitute parameters
    batch = batch.format(nodes = int(math.ceil(max(nodeslist) / pernode)))
    # List all enqueued simulations
    for i in range(len(names)):
        batch += "mpirun -np {0} -ppn {1} ./ns conf/parallel/tmp/{2}/{3}.xml 2>>output/{2}/times.txt\n".format(nodeslist[i], pernode, run_name, names[i])

    # Write new batch file
    outBatch = open("tmp/{}/batch.sh".format(run_name), 'w')
    outBatch.write(batch)
    outBatch.close()

    # Run all enqueued simulations
    subprocess.call("sbatch tmp/{}/batch.sh".format(run_name), shell=True)
