import NODE_deconvolution as nd

def read_file(path):
    f1 = open(path)
    result_file = []
    for line in f1.readlines():
        line = line.strip('\n')
        line = line.split()
        result_file.append(line)
    f1.close()
    return result_file

sc_data = read_file('E:/upload_data/SCC/filtered_data/rep2/sc_data.txt')
st_data = read_file('E:/upload_data/SCC/filtered_data/rep2/st_data.txt')
st_coordinates = read_file('E:/upload_data/SCC/filtered_data/rep2/st_coordinates.txt')
celltype = read_file('E:/upload_data/SCC/filtered_data/rep2/celltype.txt')

result_data,W_interaction = nd.get_deconvolution(
                                    st_data = st_data,
                                    sc_data = sc_data,
                                    cell_type = celltype,
                                    st_coordinate = st_coordinates,
                                    method_optimize = 1,
                                    prossecing_reserve = False,
                                    sci_max_iter = 300,
                                    file_path = 'C:/Users/HP/Desktop/test_file/test_result.txt',
                                    Number_of_iterations = 800)