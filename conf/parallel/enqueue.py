import os
import subprocess
import math

names = []
nodeslist = []
pernode = 28


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

    # Create output folder
    if not os.path.exists("../../output/{0}/{1}".format(run_name, name)):
        os.mkdir("../../output/{0}/{1}".format(run_name, name))

    # Enqueue simulation
    names.append(name)
    nodeslist.append(nProcX * nProcY * nProcZ)

def strong_scaling():
    sX = 40
    sY = 40
    sZ = 40
    lX = 1.0
    lY = 1.0
    lZ = 1.0
    enqueue_sim(cav, "s_1", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "s_2_x", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 1, 1)
    enqueue_sim(cav, "s_2_y", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 2, 1)
    enqueue_sim(cav, "s_2_z", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 2)
    enqueue_sim(cav, "s_4", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 1)
    enqueue_sim(cav, "s_8", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)
    enqueue_sim(cav, "s_16_x", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)
    # enqueue_sim(cav, "s_16_y", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 4, 2)
    # enqueue_sim(cav, "s_16_z", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 4)
    # enqueue_sim(cav, "s_32", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 4, 2)
    # enqueue_sim(cav, "s_64", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 4, 4)

def weak_scaling():
    s = 20
    lX = 1.0
    lY = 1.0
    lZ = 1.0
    enqueue_sim(cav, "w_1", "dns", "uniform", 3, lX, lY, lZ, s, s, s, 1, 1, 1)
    enqueue_sim(cav, "w_2_x", "dns", "uniform", 3, lX, lY, lZ, 2*s, s, s, 2, 1, 1)
    enqueue_sim(cav, "w_2_y", "dns", "uniform", 3, lX, lY, lZ, s, 2*s, s, 1, 2, 1)
    enqueue_sim(cav, "w_2_z", "dns", "uniform", 3, lX, lY, lZ, s, s, 2*s, 1, 1, 2)
    enqueue_sim(cav, "w_4", "dns", "uniform", 3, lX, lY, lZ, 2*s, 2*s, s, 2, 2, 1)
    enqueue_sim(cav, "w_8", "dns", "uniform", 3, lX, lY, lZ, 2*s, 2*s, 2*s, 2, 2, 2)
    enqueue_sim(cav, "w_16_x", "dns", "uniform", 3, lX, lY, lZ, 4*s, 2*s, 2*s, 4, 2, 2)

def validation():
    sX = 20
    sY = 20
    sZ = 20
    lX = 1.0
    lY = 1.0
    lZ = 1.0
    enqueue_sim(cav, "cav_2D_seq", "dns", "uniform", 2, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "cav_2D", "dns", "uniform", 2, lX, lY, lZ, sX, sY, sZ, 4, 2, 1)
    enqueue_sim(cav, "cav_3D_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "cav_3D_2", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 1, 1)
    enqueue_sim(cav, "cav_3D_8", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)
    enqueue_sim(cav, "cav_3D_16", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)
    enqueue_sim(cav, "cav_3D_64", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 4, 4)
    lY = 5.0
    enqueue_sim(cav, "cav_y5.0_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "cav_y5.0", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)
    lY = 1.0
    lZ = 0.5
    enqueue_sim(cav, "cav_z0.5_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    enqueue_sim(cav, "cav_z0.5", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 2, 2, 2)

    # sX = 20
    # sY = 20
    # sZ = 10
    # lX = 5.0
    # lY = 1.0
    # lZ = 1.0
    # enqueue_sim(bfs, "bfs_unif_seq", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    # enqueue_sim(bfs, "bfs_unif", "dns", "uniform", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)
    # enqueue_sim(bfs, "bfs_str_seq", "dns", "stretched", 3, lX, lY, lZ, sX, sY, sZ, 1, 1, 1)
    # enqueue_sim(bfs, "bfs_str", "dns", "stretched", 3, lX, lY, lZ, sX, sY, sZ, 4, 2, 2)

if __name__ == '__main__':


    # Name of this run
    run_name = "scaling_3"


    # Create folders
    if not os.path.exists("tmp/{}".format(run_name)):
        os.mkdir("tmp/{}".format(run_name))
    if not os.path.exists("../../output/{}".format(run_name)):
        os.mkdir("../../output/{}".format(run_name))


    # Enqueue simulations to run
    cav = "conf_cavity.xml"
    bfs = "conf_bfs.xml"
    strong_scaling()
    weak_scaling()
    # validation()


    # Read batch file template
    inBatch = open("batch.sh", 'r')
    batch = inBatch.read()
    inBatch.close()

    # Substitute parameters
    batch = batch.format(run_name=run_name, nodes = int(math.ceil(max(nodeslist) * 1.0 / pernode)))
    # List all enqueued simulations
    for i in range(len(names)):
        batch += "mpirun -np {0} -ppn {1} ./ns conf/parallel/tmp/{2}/{3}.xml 2>>output/{2}/times.txt\n".format(nodeslist[i], pernode, run_name, names[i])

    # Write new batch file
    outBatch = open("tmp/{}/batch.sh".format(run_name), 'w')
    outBatch.write(batch)
    outBatch.close()

    # Run all enqueued simulations
    subprocess.call("sbatch tmp/{}/batch.sh".format(run_name), shell=True)
