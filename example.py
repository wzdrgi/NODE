import NODE_deconvolution as nd
from NODE_deconvolution.NODE import get_test_data

st_data,sc_data,cell_type,st_pixel = get_test_data(10,2000,800,10)

result_data,W_interaction = nd.get_deconvolution(
                                    st_data = st_data,
                                    sc_data = sc_data,
                                    cell_type = cell_type,
                                    st_coordinate = st_pixel,
                                    method_optimize = 1,
                                    prossecing_reserve = False,
                                    file_path = '',
                                    Number_of_iterations = 500)