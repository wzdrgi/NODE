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

sc_data = read_file('E:/upload_data/developing_human_heart/filtered_data/sc_data.txt')
st_data = read_file('E:/upload_data/developing_human_heart/filtered_data/st_data_6.txt')
st_coordinates = read_file('E:/upload_data/developing_human_heart/filtered_data/st_data_coordinate_6.txt')
celltype = read_file('E:/upload_data/developing_human_heart/filtered_data/celltype_14.txt')

result_data,result_data_normalized,W_interaction = nd.get_deconvolution(
                                    st_data = st_data,
                                    sc_data = sc_data,
                                    cell_type = celltype,
                                    st_coordinate = st_coordinates,
                                    method_optimize = 1,
                                    prossecing_reserve = False,
                                    sci_max_iter = 300,
                                    file_path = 'C:/Users/HP/Desktop/test_file/test_result.txt',
                                    Number_of_iterations = 800)