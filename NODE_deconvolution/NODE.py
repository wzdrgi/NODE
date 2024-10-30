# from sko.tools import set_run_mode
from sko.DE import DE
from scipy.optimize import minimize
import random
import numpy as np
from tqdm import tqdm

#This function is used to read text files.
#Only for testing
#type_number is the maximum number of cell types present
def data_visualization(data):
    #Output the first five rows and the first five columns
    for i in range(5):
        print(data[i][0:5])
        
def get_test_data(spot_number,gene_number,cell_number,type_number):
    #Build the data first and add the table header last
    #st_data:gene * spot
    #sc_data:gene * cell
    #
    st_data = []
    for i in range(gene_number):
        list_gene = []
        for j in range(spot_number):
            gene_count = random.randrange(0,100)
            list_gene.append(str(gene_count))
        st_data.append(list_gene)


    for i in range(len(st_data)):
        st_data[i].insert(0,'gene{}'.format(i+1))
    label = []
    for k in range(len(st_data[0])):
        if k == 0:
            label.append('name')
            continue

        label.append('spot{}'.format(k))
    st_data.insert(0,label)

    sc_data = []
    for i in range(gene_number):
        list_gene = []
        for j in range(cell_number):
            gene_count = random.randrange(0,20)
            list_gene.append(str(gene_count))
        sc_data.append(list_gene)
    for i in range(len(sc_data)):
        sc_data[i].insert(0,'gene{}'.format(i+1))
    label = []
    for k in range(len(sc_data[0])):
        if k == 0:
            label.append('name')
        else:
            label.append('cell{}'.format(k))
    sc_data.insert(0,label)

    cell_type = []
    for i in range(cell_number+1):
        if i == 0:
            list_cell = []
            list_cell.append('name')
            list_cell.append('celltype')
            cell_type.append(list_cell)
            continue
        list_cell = []
        list_cell.append('cell{}'.format(i))
        ran = random.randrange(0,type_number)
        list_cell.append('type{}'.format(ran+1))
        cell_type.append(list_cell)
    
    st_pixel = [['spot','x','y']]
    for i in range(1,spot_number+1):
        list_pixel = []
        list_pixel.append('spot{}'.format(i))
        list_pixel.append(str((i//10)*1))
        list_pixel.append(str((i%10)))
        st_pixel.append(list_pixel)
    
    return st_data,sc_data,cell_type,st_pixel


def read_file(path):
    f1 = open(path)
    result_file = []
    for line in f1.readlines():
        line = line.strip('\n')
        line = line.split()
        result_file.append(line)
    f1.close()
    return result_file

#This function is used to check if sc_data files' cell names matches the cell names in celltype file.
def qc_sc_function(sc_data,cell_type):
    check_celltype = []
    check_sc_data = []

    for i in range(1,len(cell_type)):
        check_celltype.append(cell_type[i][0])
    for i in range(1,len(sc_data[0])):
        check_sc_data.append(sc_data[0][i])

    if len(check_celltype) == len(check_sc_data):
        for i in range(len(check_celltype)):
            if check_celltype[i] == check_sc_data[i]:
                continue
            else:
                print('The names of the cells do not match')
                return 0

    if len(check_celltype) != len(check_sc_data):
        print('The names of the cells do not match')
        return 0 
    
    return 1 

def qc_st_function(st_data,st_coordinate):
    check_st_coordinate = []
    check_st_data = []

    for i in range(1,len(st_coordinate)):
        check_st_coordinate.append(st_coordinate[i][0])
    for i in range(1,len(st_data[0])):
        check_st_data.append(st_data[0][i])

    if len(check_st_coordinate) == len(check_st_data):
        for i in range(len(check_st_coordinate)):
            if check_st_coordinate[i] == check_st_data[i]:
                continue
            else:
                print('The names of the spots do not match')
                return 0

    if len(check_st_coordinate) != len(check_st_data):
        print('The names of the spots do not match')
        return 0 
    
    return 1 

def qc_sc_st(sc_data,st_data):
    check_st_gene = []
    check_sc_gene = []

    for i in range(1,len(sc_data)):
        check_sc_gene.append(sc_data[i][0])
    for i in range(1,len(st_data)):
        check_st_gene.append(st_data[i][0])

    if len(check_st_gene) == len(check_sc_gene):
        for i in range(len(check_st_gene)):
            if check_st_gene[i] == check_sc_gene[i]:
                continue
            else:
                print('The names of the genes do not match')
                return 0

    if len(check_st_gene) != len(check_sc_gene):
        print('The names of the genes do not match')
        return 0 
    
    return 1 

#This function is used to read text files that requires title processing, such as sc_data file and st_data file
# def read_file_t(path):
#     f1 = open(path)
#     result_file = []
#     flag = 0
#     for line in f1.readlines():
#         if flag == 0:
#             flag = 1
#             continue
#         line = line.strip('\n')
#         line = line.split()
#         line.pop(0)
#         result_file.append(line)
#     return result_file

#This function is used to get the reference matrix to obtain the characteristics and gene expression of different cell types
def get_celltype_featrue(sc_data,cell_type):
    cell_type_dict = {}
    for i in range(1,len(cell_type)):
        cell_type_dict[cell_type[i][0]] = cell_type[i][1] 
    #In the dictionary, assign the cell type to the cell name
    type_list = []
    for k in cell_type_dict:
        if cell_type_dict[k] in type_list:
            continue
        else:
            type_list.append(cell_type_dict[k])
    #Count all cell types present
    typeshu = len(type_list)

    number_type = {}
    for number,type in enumerate(type_list):
        number_type[type] = number
    #Assigning positional numbers to cell types

    count_cell = []
    for i in range(typeshu):
        count_cell.append(0)
    for k in cell_type_dict:
        #k is cell name
        count_cell[number_type[cell_type_dict[k]]] = count_cell[number_type[cell_type_dict[k]]] + 1
    
    list_biaotou = ['name']
    for k in type_list:
        list_biaotou.append(k)
    #Add cell type information to the first line
    
    sc_mean = []
    sc_mean.append(list_biaotou)

    for i in range(1,len(sc_data)):
        list_mid = []
        for k in range(typeshu):
            list_mid.append(0)
        for j in range(1,len(sc_data[0])):
            list_mid[number_type[cell_type_dict[sc_data[0][j]]]] = list_mid[number_type[cell_type_dict[sc_data[0][j]]]] + float(sc_data[i][j])
        for m in range(typeshu):
            list_mid[m] = list_mid[m] / count_cell[m]
        list_mid.insert(0,sc_data[i][0])
        #Add gene information
        sc_mean.append(list_mid)

    return sc_mean

def get_y_matrix(st_data):
    ###spot spot1 spot2
    ###gene1 4.0  8.0
    ###gene2 0.0  2.0
    st_data = np.array(st_data)
    st_data = np.delete(st_data,0,0)
    st_data = np.delete(st_data,0,1)
    st_data = st_data.T
    st_data = np.array(st_data,dtype=np.float32)
    #spot * gene
    return st_data

def get_beta_matrix(sc_mean):
    sc_mean = np.array(sc_mean)
    sc_mean = np.delete(sc_mean,0,0)
    sc_mean = np.delete(sc_mean,0,1)
    sc_mean = sc_mean.T
    sc_mean = np.array(sc_mean,dtype=np.float32)
    #type * gene
    return sc_mean

def get_d_matrix(st_coordinate,neighbors_number):
    ### spot    x       y
    ### spot1   8       20
    ### spot2   10      12
    ### 会改变st_coordinate
    # st_coordinate.pop(0)
    spot_number = len(st_coordinate) - 1
    d = np.zeros((spot_number,spot_number))
    for k in range(spot_number):
        for j in range(spot_number):
            d[k][j] = np.sqrt((float(st_coordinate[k+1][1])-float(st_coordinate[j+1][1]))**2 + (float(st_coordinate[k+1][2])-float(st_coordinate[j+1][2]))**2)
    #Also calculated the distance from itself to itself
    #Take the first four SPOTS and, since it calculates itself, take 5
    list_lim = []
    for i in d:
        k = np.sort(i)
        list_lim.append(k[neighbors_number]) 

    for i in range(spot_number):
        for j in range(spot_number):
            if d[i][j] > list_lim[i]:
                d[i][j] = 0
            if d[i][j] <= list_lim[i] and d[i][j] != 0:
                d[i][j] = 1
    d = np.array(d,dtype=np.float32)
    return d

def Remove_low_quality_ST(ST_data,st_coordinare):
    #首先获得ST_data的大小
    #要求是代表头的方阵
    a = len(ST_data)
    b = len(ST_data[0])
    filterd_list = []
    for i in range(1,b):
        count = 0
        for j in range(1,a):
            count = count + float(ST_data[j][i])

        if count < 100:
            filterd_list.append(i)

        else:
            continue
    
    filterd_list.sort(reverse=True)
    print(filterd_list)
    for i in range(a):
        for j in filterd_list:
            ST_data[i].pop(j)
    
    st_coordinare_new = []
    for i in range(len(st_coordinare)):
        if i in filterd_list:
            continue
        else:
            st_coordinare_new.append(st_coordinare[i])
    
    return ST_data,st_coordinare_new

#sc_data,st_data,celltype,st_coordinate are lists
#In the main function of the deconvolution, only three data are needed, which are sc_data,st_data,sc_mean.
def get_STinteractions(expre_matrix, Beta_matrix, deconvolution_x):
    # Y = W * Y + X * B
    expre_matrix = np.array(expre_matrix, dtype = np.float32)
    # spot * gene
    Beta_matrix = np.array(Beta_matrix, dtype = np.float32)
    # type * gene
    deconvolution_x = np.array(deconvolution_x, dtype = np.float32)
    # spot * type
    cell_expre = np.dot(deconvolution_x, Beta_matrix)
    interaction_expre = expre_matrix - cell_expre
    expre_matrix_inverse = np.linalg.pinv(expre_matrix)
    W_communication = np.dot(interaction_expre,expre_matrix_inverse)
    # print(W_communication.shape)
    return W_communication

def interaction_communication(W_interaction_deconvolution,D_matrix,absolute_intensity):
    W_interaction_deconvolution = np.array(W_interaction_deconvolution, dtype = np.float32)
    W_interaction_deconvolution = W_interaction_deconvolution.T
    # Row spot sends out signals, column spot receives signals, 
    # and positive and negative represent the up- and down-regulation of the signal on the gene expression of the target spot
    if absolute_intensity:
        #Enable absolute intensity, then do not show positive or negative in the
        W_interaction_deconvolution[W_interaction_deconvolution < 0] = -W_interaction_deconvolution[W_interaction_deconvolution < 0]
    W_interaction_deconvolution = np.multiply(W_interaction_deconvolution,D_matrix)
    return W_interaction_deconvolution

def add_lablenames(W_interaction_communication,st_coordinate):
    W_interaction_scipy_communication = W_interaction_communication.tolist()
    spot_name = []
    for i in range(1,len(st_coordinate)):
        spot_name.append(st_coordinate[i][0])
        W_interaction_scipy_communication[i-1].insert(0,st_coordinate[i][0])
    spot_name.insert(0,'spot')
    W_interaction_scipy_communication.insert(0,spot_name)
    return W_interaction_scipy_communication

def get_deconvolution(sc_data,
                      st_data, 
                      cell_type, 
                      st_coordinate, 
                      neighbors_number = 4,
                      #Number of neighbors considered
                      method_optimize = 2,
                      prossecing_reserve = False,
                      #Recommend saving while deconvolution
                      #Whether or not it is saved in the process
                      file_path = '',
                      #Save Path
                      sci_max_iter = 150,
                      max_cell = 8,
                      max_iter = 800,
                      #Maximum number of iterations
                      Number_of_iterations = 700,
                      max_interaction = 1,
                      #Actual number of iterations
                      spot_start = 0,
                      #If you are resuming after an interruption, you need to update the save path, otherwise the original file will be overwritten.
                      absolute_intensity = False,
                      #The output form of the communication matrix, the rows and columns of the communication matrix represent the communication line of defense, 
                      #and the positive and negative matrix elements represent the upward and downward adjustments; 
                      #if this value is set to True, then the upward and downward adjustments are not displayed, 
                      #and only the communication strength is outputted.
):
    #Creating chokepoints and deciding whether to enter the deconvolution program or not
    #STEP1: QC
    if prossecing_reserve != True and prossecing_reserve != False:
        print("Save process boolean set incorrectly. The save process boolean is set incorrectly and can only be True or False.")
        return 0

    if isinstance(neighbors_number,int) != True:
        print("neighbors_number can only be of type int.")
        return 0 
    
    if method_optimize != 1 and method_optimize != 2:
        print("Please check the method_optimize, which can only be 1 or 2")
        return 0
    
    if isinstance(sci_max_iter,int) != True:
        print("sci_max_iter can only be of type int.")
        return 0 
    
    if isinstance(max_cell,int) != True:
        print("max_cell can only be of type int.")
        return 0 

    if isinstance(max_iter,int) != True:
        print("max_iter can only be of type int")
        return 0
    
    if isinstance(Number_of_iterations,int) != True:
        print("Number_of_iterations can only be of type int")
        return 0
    
    if max_interaction > 1 or max_interaction < 0:
        print("Please set max_interaction between 0 and 1")
        return 0
    
    if isinstance(spot_start,int) != True:
        print("spot_start can only be of type int")
        return 0

    if spot_start > len(st_data) - 2:
        print("The spot start position cannot exceed the total number of spots")
        return 0

    if prossecing_reserve == True:
        if file_path == '':
            print("Please set the save path during the work process")
            return 0
    
    if neighbors_number + 1 >= len(st_coordinate):
        print("Selection of neighbor spots exceeds all the number of spots, Please don't choose too many neighboring spots")
        return 0
        
    if prossecing_reserve == True:
        fun_flag = 1 
        if qc_sc_function(sc_data = sc_data,cell_type = cell_type) == 0:
            fun_flag = 0

        if qc_st_function(st_data=st_data,st_coordinate=st_coordinate) == 0:
            fun_flag = 0

        if qc_sc_st(sc_data=sc_data,st_data=st_data) == 0:
            fun_flag = 0

        if fun_flag == 0:
            print("please check the information of input data")
            return 0
        #Processes the data so that it can be passed as a parameter to the
        #STEP2: creat sc_mean

        sc_mean = get_celltype_featrue(sc_data=sc_data,cell_type=cell_type)
        sc_type_lable = sc_mean[0]
        st_lable = []

        for i in range(1,len(st_data[0])):
            st_lable.append(st_data[0][i])
        #STEP3: creat three parameters
        y_matrix = get_y_matrix(st_data=st_data)
        #spot * gene
        beta_matrix = get_beta_matrix(sc_mean=sc_mean)
        #type * gene
        d_matrix = get_d_matrix(st_coordinate=st_coordinate,neighbors_number=neighbors_number)
        #STEP4: deconvolution
        #First set the solution mode

        solution_mode = {1:'sci',2:'de_model'}

        if solution_mode[method_optimize] == 'sci':
            ##use scipy
            print('use scipy to deconvolution')
            print('Deconvolution Start!') 
            result_cell = []
            result_cell.append(sc_type_lable)

            type_number = len(beta_matrix)
            spot_number = len(y_matrix)

            # print(y_matrix.shape)
            # print(beta_matrix.shape)
            # print(d_matrix.shape)
            # print(type_number)
            # print(spot_number)

            for spot in tqdm(range(spot_start,spot_number)):
                y_matrix_test = y_matrix[spot,:]
                W0 = np.zeros((1,spot_number))
                def fun_x_0(X):
                    return  np.linalg.norm((np.multiply(W0,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #W is a constant that optimizes X
                cons1 = []
                for i in range(type_number):
                    cons1.append({'type': 'ineq', 'fun': lambda X: X[i]})

                X0 = np.zeros(type_number)
                res = minimize(fun_x_0, X0, method='SLSQP' , constraints=cons1, options={'maxiter': sci_max_iter,'disp': False})
                x_mid = res.x
                def fun_w_0(W):
                    X = x_mid
                    return np.linalg.norm((np.multiply(W,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #X is a constant optimizing W
                cons = []
                for i in range(spot_number):
                    cons.append({'type': 'ineq', 'fun': lambda W: -W[i ]+ max_interaction})
                    cons.append({'type': 'ineq', 'fun': lambda W: W[i] + max_interaction})
                W0 = np.zeros(spot_number)
                res = minimize(fun_w_0, W0, method='SLSQP' , constraints=cons, options={'maxiter': sci_max_iter,'disp': False})
                W_mid = res.x

                def fun_x(X):
                    W = W_mid
                    return  np.linalg.norm( np.multiply(d_matrix[spot],W) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #W is a constant that optimizes X
                def fun_w(W):
                    X = x_mid
                    return np.linalg.norm(np.multiply(d_matrix[spot],W) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #X is a constant optimizing W

                for i in range(1):
                    X0 = x_mid
                    res = minimize(fun_x, X0, method='SLSQP' , constraints=cons1, options={'maxiter': sci_max_iter,'disp': False})
                    x_mid = res.x
                
                save_list = []
                save_list.append(st_lable[spot])
                for k in x_mid:
                    if float(k) < 0:
                        save_list.append(0)
                        continue
                    if float(k) > max_cell:
                        save_list.append(max_cell)
                        continue
                    
                    save_list.append(float(k))
                #guaranteed non-negative
                result_cell.append(save_list)
                np.savetxt(file_path,result_cell,fmt='%s')
                print('Saving during the process{}/{}'.format(spot+1,spot_number))
            
            result_cell = np.array(result_cell)
            result_x_mid = result_cell[1:,1:]

            W_ini = get_STinteractions(expre_matrix = y_matrix,Beta_matrix = beta_matrix,deconvolution_x = result_x_mid)
            W_interaction = interaction_communication(W_ini,d_matrix,absolute_intensity)

            result_cell_normalized = result_cell.copy()
            result_cell_normalized_array = result_cell_normalized[1:,1:]
            result_cell_normalized_array = result_cell_normalized_array.astype(float)
            row_sums = result_cell_normalized_array.sum(axis=1, keepdims=True) 
            row_sums[row_sums == 0] = 1
            normalized_data = result_cell_normalized_array / row_sums
            result_cell_normalized[1:, 1:] = normalized_data

            result_cell = result_cell.tolist()
            result_cell_normalized = result_cell_normalized.tolist()

            W_interaction = add_lablenames(W_interaction,st_coordinate=st_coordinate)

            print("Deconvolution Finish!")
            return result_cell,result_cell_normalized,W_interaction

        

        if solution_mode[method_optimize] == 'de_model':
            result = []
            result.append(sc_type_lable)
            print('use scikit-opt to deconvolution')
            print('Deconvolution Start') 
            type_number = len(beta_matrix)
            spot_number = len(y_matrix)

            def obj_func(X):
                return  np.linalg.norm( (np.multiply(W0,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
            #First optimize X
            # mode = 'multithreading'
            #multiprocess acceleration
            # set_run_mode(obj_func, mode)

            '''
            Upper and lower bounds for initial x
            '''
            ini_x_lim = []
            ini_x_up = []
            ini_w_lim = []
            ini_w_up = []

            for i in range(spot_number):
                ini_w_lim.append(-1 * max_interaction)
            for i in range(spot_number):
                ini_w_up.append(1 * max_interaction)
            for i in range(type_number):
                ini_x_lim.append(0)
            for i in range(type_number):
                ini_x_up.append(max_cell)
            


            for spot in tqdm(range(spot_start,spot_number)):

                y_matrix_test = y_matrix[spot,:]
                W0 = np.zeros((1,spot_number))

                constraint_ueq = []
                #restrictive condition,The constraints come from the cell types, so the conditions are the same as the number of cell types
                constraint_ueq_w = []
                for i in range(spot_number):
                    constraint_ueq_w.append(lambda W: -W[i])
                for i in range(type_number):
                    constraint_ueq.append(lambda x: -x[i])


                de = DE(func = obj_func, n_dim= type_number, size_pop=50, max_iter = max_iter, lb=ini_x_lim, ub=ini_x_up,constraint_ueq=constraint_ueq)
                best_x, best_y = de.run(Number_of_iterations)
                #best_x : parameter resolution
                #best_y : objective function value
                X_mid = best_x
                array_1 = best_x


                def fun_w_0(W):
                    X = best_x
                    return np.linalg.norm((np.multiply(W,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                # set_run_mode(fun_w_0, mode)
                de = DE(func = fun_w_0, n_dim= spot_number, size_pop=50, max_iter = max_iter, lb = ini_w_lim, ub=ini_w_up, constraint_ueq=constraint_ueq_w)
                best_x, best_y = de.run(200)

                def obj_func_1(X):
                    W0 = best_x
                    return  np.linalg.norm( (np.multiply(W0,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                # set_run_mode(obj_func_1, mode)
                X_mid = X_mid * 0.8
                de = DE(func = obj_func_1, n_dim = type_number, size_pop=50, max_iter = max_iter, lb= X_mid, ub=ini_x_up,constraint_ueq=constraint_ueq)
                best_x, best_y = de.run(Number_of_iterations)
                list_mid = []
                list_mid.append(st_lable[spot])
                for i in best_x:
                    list_mid.append(i)
                result.append(list_mid)

                
                np.savetxt(file_path,result,fmt='%s')
                print('Saving during the process{}/{}'.format(spot+1,spot_number))

            result = np.array(result)
            result_x_mid = result[1:,1:]

            W_ini = get_STinteractions(expre_matrix = y_matrix,Beta_matrix = beta_matrix,deconvolution_x = result_x_mid)
            W_interaction = interaction_communication(W_ini,d_matrix,absolute_intensity)
            
            result_normalized = result.copy()
            result_normalized_array = result_normalized[1:,1:]
            result_normalized_array = result_normalized_array.astype(float)
            row_sums = result_normalized_array.sum(axis=1,keepdims=True)
            row_sums[row_sums == 0] = 1
            normalized_data = result_normalized_array / row_sums
            result_normalized[1:, 1:] = normalized_data

            result = result.tolist()
            result_normalized = result_normalized.tolist()
            
            W_interaction = add_lablenames(W_interaction,st_coordinate=st_coordinate)
            
            print("Deconvolution Finish!")
            return result,result_normalized,W_interaction
    else:

        fun_flag = 1 
        if qc_sc_function(sc_data = sc_data,cell_type = cell_type) == 0:
            fun_flag = 0
        if qc_st_function(st_data=st_data,st_coordinate=st_coordinate) == 0:
            fun_flag = 0
        if qc_sc_st(sc_data=sc_data,st_data=st_data) == 0:
            fun_flag = 0
        if fun_flag == 0:
            print("please check the information of input data")
            return 0
        #Processes the data so that it can be passed as a parameter to the
        #STEP2: creat sc_mean
        sc_mean = get_celltype_featrue(sc_data=sc_data,cell_type=cell_type)
        sc_type_lable = sc_mean[0]
        st_lable = []
        for i in range(1,len(st_data[0])):
            st_lable.append(st_data[0][i])
        #STEP3: creat three parameters
        y_matrix = get_y_matrix(st_data=st_data)
        #spot * gene
        beta_matrix = get_beta_matrix(sc_mean=sc_mean)
        d_matrix = get_d_matrix(st_coordinate=st_coordinate,neighbors_number=neighbors_number)
        #STEP4: deconvolution
        #First set the solution mode
        solution_mode = {1:'sci',2:'de_model'}

        if solution_mode[method_optimize] == 'sci':
            ##use scipy
            print('use scipy to deconvolution')
            print('Deconvolution Start') 
            result_cell = []
            result_cell.append(sc_type_lable)
            type_number = len(beta_matrix)
            spot_number = len(y_matrix)

            # print(y_matrix.shape)
            # print(beta_matrix.shape)
            # print(d_matrix.shape)
            # print(type_number)
            # print(spot_number)

            for spot in tqdm(range(spot_start,spot_number)):
                y_matrix_test = y_matrix[spot,:]
                W0 = np.zeros((1,spot_number))
                def fun_x_0(X):
                    return  np.linalg.norm((np.multiply(W0,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #W is a constant that optimizes X
                cons1 = []
                for i in range(type_number):
                    cons1.append({'type': 'ineq', 'fun': lambda X: X[i]})

                X0 = np.zeros(type_number)
                res = minimize(fun_x_0, X0, method='SLSQP' , constraints=cons1, options={'maxiter': sci_max_iter,'disp': False})
                x_mid = res.x
                def fun_w_0(W):
                    X = x_mid
                    return np.linalg.norm((np.multiply(W,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #X is a constant optimizing W
                cons = []
                for i in range(spot_number):
                    cons.append({'type': 'ineq', 'fun': lambda W: -W[i ]+ max_interaction})
                    cons.append({'type': 'ineq', 'fun': lambda W: W[i] + max_interaction})
                W0 = np.zeros(spot_number)
                res = minimize(fun_w_0, W0, method='SLSQP' , constraints=cons, options={'maxiter': sci_max_iter,'disp': False})
                W_mid = res.x

                def fun_x(X):
                    W = W_mid
                    return  np.linalg.norm( np.multiply(d_matrix[spot],W) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #W is a constant that optimizes X
                def fun_w(W):
                    X = x_mid
                    return np.linalg.norm(np.multiply(d_matrix[spot],W) @ y_matrix + X @ beta_matrix - y_matrix_test)
                #X is a constant optimizing W

                for i in range(1):
                    X0 = x_mid
                    res = minimize(fun_x, X0, method='SLSQP' , constraints=cons1, options={'maxiter': sci_max_iter,'disp': False})
                    x_mid = res.x
                
                save_list = []
                save_list.append(st_lable[spot])
                for k in x_mid:
                    if float(k) < 0:
                        save_list.append(0)
                        continue
                    if float(k) > max_cell:
                        save_list.append(max_cell)
                        continue

                    save_list.append(float(k))
                #guaranteed non-negative
                result_cell.append(save_list)

            result_cell = np.array(result_cell)
            result_x_mid = result_cell[1:,1:]

            W_ini = get_STinteractions(expre_matrix = y_matrix,Beta_matrix = beta_matrix,deconvolution_x = result_x_mid)
            W_interaction = interaction_communication(W_ini,d_matrix,absolute_intensity)

            result_cell_normalized = result_cell.copy()
            result_cell_normalized_array = result_cell_normalized[1:,1:]
            result_cell_normalized_array = result_cell_normalized_array.astype(float)
            row_sums = result_cell_normalized_array.sum(axis=1, keepdims=True) 
            row_sums[row_sums == 0] = 1
            normalized_data = result_cell_normalized_array / row_sums
            result_cell_normalized[1:, 1:] = normalized_data

            result_cell = result_cell.tolist()
            result_cell_normalized = result_cell_normalized.tolist()

            W_interaction = add_lablenames(W_interaction,st_coordinate=st_coordinate)

            print("Deconvolution Finish!")
            return result_cell,result_cell_normalized,W_interaction

        

        if solution_mode[method_optimize] == 'de_model':
            result = []
            result.append(sc_type_lable)
            print('use scikit-opt to deconvolution')
            print('Deconvolution Start') 
            type_number = len(beta_matrix)
            spot_number = len(y_matrix)

            def obj_func(X):
                return  np.linalg.norm( (np.multiply(W0,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
            #First optimize X
            # mode = 'multithreading'
            #multiprocess acceleration
            # set_run_mode(obj_func, mode)

            '''
            Upper and lower bounds for initial x
            '''
            ini_x_lim = []
            ini_x_up = []
            ini_w_lim = []
            ini_w_up = []


            for i in range(spot_number):
                ini_w_lim.append(-1 * max_interaction)
            for i in range(spot_number):
                ini_w_up.append(1 * max_interaction)
            for i in range(type_number):
                ini_x_lim.append(0)
            for i in range(type_number):
                ini_x_up.append(max_cell)
            


            for spot in tqdm(range(spot_start,spot_number)):

                y_matrix_test = y_matrix[spot,:]
                W0 = np.zeros((1,spot_number))

                constraint_ueq = []
                #restrictive condition,The constraints come from the cell types, so the conditions are the same as the number of cell types
                constraint_ueq_w = []
                for i in range(spot_number):
                    constraint_ueq_w.append(lambda W: -W[i])
                for i in range(type_number):
                    constraint_ueq.append(lambda x: -x[i])


                de = DE(func = obj_func, n_dim= type_number, size_pop=50, max_iter = max_iter, lb=ini_x_lim, ub=ini_x_up,constraint_ueq=constraint_ueq)
                best_x, best_y = de.run(Number_of_iterations)
                #best_x : parameter resolution
                #best_y : objective function value
                X_mid = best_x
                array_1 = best_x


                def fun_w_0(W):
                    X = best_x
                    return np.linalg.norm((np.multiply(W,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                # set_run_mode(fun_w_0, mode)
                de = DE(func = fun_w_0, n_dim= spot_number, size_pop=50, max_iter = max_iter, lb = ini_w_lim, ub=ini_w_up, constraint_ueq=constraint_ueq_w)
                best_x, best_y = de.run(200)

                def obj_func_1(X):
                    W0 = best_x
                    return  np.linalg.norm( (np.multiply(W0,d_matrix[spot])) @ y_matrix + X @ beta_matrix - y_matrix_test)
                # set_run_mode(obj_func_1, mode)
                X_mid = X_mid * 0.8
                de = DE(func = obj_func_1, n_dim = type_number, size_pop=50, max_iter = max_iter, lb= X_mid, ub=ini_x_up,constraint_ueq=constraint_ueq)
                best_x, best_y = de.run(Number_of_iterations)
                list_mid = []
                list_mid.append(st_lable[spot])
                for i in best_x:
                    list_mid.append(i)
                result.append(list_mid)
            
            result = np.array(result)
            result_x_mid = result[1:,1:]

            W_ini = get_STinteractions(expre_matrix = y_matrix,Beta_matrix = beta_matrix,deconvolution_x = result_x_mid)
            W_interaction = interaction_communication(W_ini,d_matrix,absolute_intensity)
            
            result_normalized = result.copy()
            result_normalized_array = result_normalized[1:,1:]
            result_normalized_array = result_normalized_array.astype(float)
            row_sums = result_normalized_array.sum(axis=1,keepdims=True)
            row_sums[row_sums == 0] = 1
            normalized_data = result_normalized_array / row_sums
            result_normalized[1:, 1:] = normalized_data

            result = result.tolist()
            result_normalized = result_normalized.tolist()
            
            W_interaction = add_lablenames(W_interaction,st_coordinate=st_coordinate)
            
            print("Deconvolution Finish!")
            return result,result_normalized,W_interaction
