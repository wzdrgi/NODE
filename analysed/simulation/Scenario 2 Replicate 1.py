import numpy as np
from sko.DE import DE
from tqdm import tqdm

def read_file(path):
    f1 = open(path)
    sc_data = []
    for line in f1.readlines():
        line = line.strip('\n')
        line = line.split()
        sc_data.append(line)
    f1.close()
    return sc_data

y_matrix_ini = read_file('E:/upload_data/simulation/Scenario 2 Replicate 1/st_data.txt')
y_matrix_ini.pop(0)
for i in range(len(y_matrix_ini)):
    y_matrix_ini[i].pop(0)
y_matrix = np.array(y_matrix_ini,dtype=np.float32)

d_matrix_ini = read_file('E:/upload_data/simulation/Scenario 2 Replicate 1/st_coordinates.txt')
d_matrix_ini.pop(0)
for i in range(len(d_matrix_ini)):
    d_matrix_ini[i].pop(0)
spot_number = len(d_matrix_ini)

d = np.zeros((spot_number,spot_number))
for k in range(spot_number):
    for j in range(spot_number):
        d[k][j] = np.sqrt((float(d_matrix_ini[k][0])-float(d_matrix_ini[j][0]))**2 + (float(d_matrix_ini[k][1])-float(d_matrix_ini[j][1]))**2)

list_lim = []
for i in d:
    k = np.sort(i)
    list_lim.append(k[5]) 

for i in range(spot_number):
    for j in range(spot_number):
        if d[i][j] >= list_lim[i]:
            d[i][j] = 0
        if d[i][j] < list_lim[i] and d[i][j] != 0:
            d[i][j] = 1
d = np.array(d)

beta_matrix = read_file('E:/upload_data/simulation/Scenario 2 Replicate 1/sc_mean.txt')
beta_matrix.pop(0)
for i  in range(len(beta_matrix)):
    beta_matrix[i].pop(0)
beta_matrix = np.array(beta_matrix,dtype=np.float32)

beta_matrix = beta_matrix.T
y_matrix = y_matrix.T

confusion_matrix = []
for i in y_matrix:
    i = i.tolist()
    confusion_matrix.append(i)
for i in beta_matrix:
    i = i.tolist()
    confusion_matrix.append(i)

confusion_matrix = np.array(confusion_matrix,dtype=np.float32)

spot_number = len(y_matrix)
type_number = len(beta_matrix)
decision_variables_number = spot_number + type_number

result = []
for spot in tqdm(range(spot_number)):

    y_matrix_test = y_matrix[spot,:]

    new_d = []
    for i in d[spot]:
        new_d.append(i)
    for i in range(type_number):
        new_d.append(1)
    new_d = np.array(new_d,dtype=np.float32)

    def obj_func(decision_variables):
        return  np.linalg.norm((np.multiply(decision_variables,new_d)) @ confusion_matrix - y_matrix)

    qianzhi_lb = []
    for i in range(spot_number):
        qianzhi_lb.append(0)
    for i in range(type_number):
        qianzhi_lb.append(0)

    houzhi_ub = []
    for i in range(spot_number):
        if new_d[i] == 0:
            houzhi_ub.append(0)
        else:
            houzhi_ub.append(1)
    for i in range(type_number):
        houzhi_ub.append(8)

    constraint_ueq_w = []
    for i in range(spot_number):
        constraint_ueq_w.append(lambda W: -W[i])
    for i in range(type_number):
        constraint_ueq_w.append(lambda W: -W[i])

    de = DE(func = obj_func, n_dim = spot_number + type_number, 
    size_pop = 60, max_iter = 1000, 
    lb = qianzhi_lb, ub=houzhi_ub, constraint_ueq=constraint_ueq_w)
    best_x, best_y = de.run(1000)
    result.append(best_x)

result = np.array(result)
result_W = result[:,0:len(result)-type_number]
result_X = result[:,len(result)-type_number:]

np.savetxt('E:/upload_data/simulation/result/Scenario 2 Replicate 1.txt',result_X,fmt='%s')
np.savetxt('E:/upload_data/simulation/result/Scenario 2 Replicate 1 CCC.txt',result_W,fmt='%s')
