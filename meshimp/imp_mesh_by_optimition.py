import argparse
import torch
from torch.autograd import Variable # torch 中 Variable 模块
import torch.optim as optim
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import meshplex
from scipy.spatial import Delaunay
import numpy as np
from sklearn.externals import joblib
import time
import meshio
import meshplex
import trimesh


def laplacian(ring):
    newpoints = np.mean(ring, axis=0)
    return newpoints

def calc_area(p1, p2, p3):
    p4 = p2 - p1
    p5 = p3 - p1

    v = torch.abs(p4[0]*p5[1]-p4[1]*p5[0])

    return v

def calc_len(p1, p2, p3):
    p4 = p2 - p1
    p5 = p3 - p2
    p6 = p1 - p3

    v1 = p4.dot(p4)
    v2 = p5.dot(p5)
    v3 = p6.dot(p6)

    v = torch.sqrt(v1) + torch.sqrt(v2) + torch.sqrt(v3)

    return v

def calc_len2(p1, p2, p3):
    p4 = p2 - p1
    p5 = p3 - p2
    p6 = p1 - p3

    v1 = p4.dot(p4)
    v2 = p5.dot(p5)
    v3 = p6.dot(p6)

    v = v1 + v2 + v3

    return v

def compute_ele_energy(cells_index,tfpoints,cells):
    energy = 0
    # 计算每个单元的质量
    for e in cells_index:
        if e == -1:
            continue
        else:
            C0 = tfpoints[cells[e,0]]
            C1 = tfpoints[cells[e,1]]
            C2 = tfpoints[cells[e,2]]
            len = calc_len2(C0, C1, C2)
            area = calc_area(C0, C1, C2)
            energy = energy + len / area
    return energy

def output(ring, f):
    s = ''
    for i, num in enumerate(ring):
        s += '{} '.format(num)
    s += '\n'
    f.write(s)

def run(input_file, output_file, opt_epoch, lr, coltrol):

    #读入网格
    mesh = trimesh.load(input_file)   #这个网格结构可以方便的寻找点面之间的邻接关系
    points = mesh.vertices    #网格的点的坐标
    faces = mesh.faces   #网格的三角片
    meshbac = meshplex.MeshTri(points,faces) #可以用来查找哪些点是边界点
    bpindex = meshbac.is_boundary_node   #这个是bool数组
    vvlist = mesh.vertex_neighbors  #点的邻接点关系
    vfarray = mesh.vertex_faces

    #meshbac.show()   #展示优化前网格

    #创建变量，将所有节点都设置成变量
    variapoints = []
    for i in range(0, points.shape[0]):
        torch_point = torch.from_numpy(points[i])
        tfpoint = Variable(torch_point, requires_grad=True)
        variapoints.append(tfpoint)

    #创建输出文件

    f3 = open("./data/dataset3.txt", 'a')
    f4 = open("./data/dataset4.txt", 'a')
    f5 = open("./data/dataset5.txt", 'a')
    f6 = open("./data/dataset6.txt", 'a')
    f7 = open("./data/dataset7.txt", 'a')
    f8 = open("./data/dataset8.txt", 'a')
    f9 = open("./data/dataset9.txt", 'a')
    gt3 = open("./data/gtset3.txt", 'a')
    gt4 = open("./data/gtset4.txt", 'a')
    gt5 = open("./data/gtset5.txt", 'a')
    gt6 = open("./data/gtset6.txt", 'a')
    gt7 = open("./data/gtset7.txt", 'a')
    gt8 = open("./data/gtset8.txt", 'a')
    gt9 = open("./data/gtset9.txt", 'a')
    #print('%d,%d' % (len(matches_a), len(matches_b)), file=f)

    #优化网格
    start = time.time()
    for t in range(opt_epoch):
        print("this is epoch ",t)
        for i in range(0, points.shape[0]):
            if bpindex[i] == True:
                continue
            else:
                optimizer = optim.Adam([variapoints[i]], lr=lr)
                energy = 0
                for input in range(1, 5000):
                    optimizer.zero_grad()
                    energy1 = compute_ele_energy(vfarray[i], variapoints, faces)
                    if abs(energy1 - energy) < coltrol:
                        #print(input)
                        break
                    energy = energy1
                    #print(energy)
                    energy.backward()
                    optimizer.step()
            #保存优化结果
            '''
            if len(vvlist[i])==3:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f3)
                output(points[i], gt3)
            elif len(vvlist[i])==4:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f4)
                output(points[i], gt4)
            elif len(vvlist[i])==5:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f5)
                output(points[i], gt5)
            elif len(vvlist[i])==6:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f6)
                output(points[i], gt6)
            elif len(vvlist[i])==7:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f7)
                output(points[i], gt7)
            elif len(vvlist[i])==8:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f8)
                output(points[i], gt8)
            elif len(vvlist[i])==9:
                output(np.asarray(points[vvlist[i]]).reshape(-1),f9)
                output(points[i], gt9)
        #保存网格
        #meshbac.save("./opt2000-%d.stl"%t)
        '''
    meshbac.save(output_file)
    #关闭输出文件
    f3.close()
    f4.close()
    f5.close()
    f6.close()
    f7.close()
    f8.close()
    f9.close()
    gt3.close()
    gt4.close()
    gt5.close()
    gt6.close()
    gt7.close()
    gt8.close()
    gt9.close()
    end = time.time()
    print("this time ", t, " is ",end - start)
    meshbac.show()    #展示优化后网格

def parse_args():
    """ Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="optimize mesh by optimition")
    parser.add_argument(
        "--input_file", help="The input file.", default=None,
        required=True)
    parser.add_argument(
        "--output_file", help="The output file.", default=None,
        required=True)
    parser.add_argument(
        "--opt_epoch", help="the epoch of smooth.", type=int, default=10)
    parser.add_argument(
        "--lr", help="the learning rate.", type=float, default=0.001)
    parser.add_argument(
        "--cl", help="the control number.", type=float, default=0.000000000000000001)
    return parser.parse_args()
'''
if __name__ == "__main__":
    args = parse_args()
    run(args.input_file, args.output_file, args.opt_epoch, args.lr, args.cl)
'''
run("./sresult3.stl","./optsresult3.stl",10,0.1,0.00001)