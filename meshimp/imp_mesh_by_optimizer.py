import torch
import argparse
import numpy as np
from sklearn.externals import joblib
import time
import meshplex
import trimesh

def run(input_file, output_file, opt_epoch):
    #创建优化器
    opt3 = torch.load("./modeldir/opt3.pkl")
    opt4 = torch.load("./modeldir/opt4.pkl")
    opt5 = torch.load("./modeldir/opt5.pkl")
    opt6 = torch.load("./modeldir/opt6.pkl")
    opt7 = torch.load("./modeldir/opt7.pkl")
    opt8 = torch.load("./modeldir/opt8.pkl")
    opt9 = torch.load("./modeldir/opt9.pkl")


    def improveRing(ring):
        if ring.shape[0] == 6:
            opt = opt3
        if ring.shape[0] == 8:
            opt = opt4
        if ring.shape[0] == 10:
            opt = opt5
        if ring.shape[0] == 12:
            opt = opt6
        if ring.shape[0] == 14:
            opt = opt7
        if ring.shape[0] == 16:
            opt = opt8
        if ring.shape[0] == 18:
            opt = opt9
        ring = torch.from_numpy(ring).float()
        opoint = opt(ring).data.numpy()
        return opoint

    def laplacian(ring):
        newpoints = np.mean(ring, axis=0)
        return newpoints

    # 读入网格
    mesh = trimesh.load(input_file)  # 这个网格结构可以方便的寻找点面之间的邻接关系
    points = mesh.vertices  # 网格的点的坐标
    faces = mesh.faces  # 网格的三角片
    meshbac = meshplex.MeshTri(points, faces)  # 可以用来查找哪些点是边界点
    bpindex = meshbac.is_boundary_node  # 这个是bool数组
    vvlist = mesh.vertex_neighbors  # 点的邻接关系
    print("hello1")
    meshbac.show()  # 展示优化前网格


    # 优化网格
    start = time.time()
    for t in range(opt_epoch):
        for i in range(0, points.shape[0]):
            if bpindex[i] == True:
                continue
            else:
                vvpoints = points[vvlist[i]]
                if vvpoints.shape[0] < 3 or vvpoints.shape[0] > 9:
                    points[i]=laplacian(vvpoints)
                else:
                    vvpoints = vvpoints.flatten()
                    l = []
                    for m in range(0, vvpoints.shape[0]):
                        if (m + 1) % 3 != 0:
                            l.append(m)
                    vvpoints = vvpoints[l]

                    ringx = vvpoints[0:vvpoints.shape[0]:2]
                    ringy = vvpoints[1:vvpoints.shape[0]:2]

                    ringxmin = ringx.min()
                    ringxmax = ringx.max()
                    ringymin = ringy.min()
                    ringymax = ringy.max()

                    if ringxmax - ringxmin > ringymax-ringymin:
                        len = ringxmax - ringxmin
                    else:
                        len = ringymax - ringymin

                    vvpoints[0:vvpoints.shape[0]:2] = (ringx - ringxmin) / (len)
                    vvpoints[1:vvpoints.shape[0]:2] = (ringy - ringymin) / (len)

                    points[i, 0:2] = improveRing(vvpoints)*len+[ringxmin, ringymin]

        end = time.time()
        print("this time ", t, " is ",end - start)
        #meshbac.save("./myoptbadtri-%d.stl"%t)
    # 保存网格
    meshbac.save(output_file)
    print("hello2")
    meshbac.show()  # 展示优化后网格

def parse_args():
    """ Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="optimize mesh by optimizer")
    parser.add_argument(
        "--input_file", help="The input file.", default=None,
        required=True)
    parser.add_argument(
        "--output_file", help="The output file.", default=None,
        required=True)
    parser.add_argument(
        "--opt_epoch", help="the epoch of smooth.", type=int, default=10)
    return parser.parse_args()


'''
if __name__ == "__main__":
    args = parse_args()
    run(args.input_file, args.output_file, args.opt_epoch)
'''
run("./sresult3.stl","./myoptsresult3.stl",10)
#python imp_mesh_by_optimizer.py --input_file=./data/badtri.stl --output_file=./optbadtri.stl --opt_epoch=50