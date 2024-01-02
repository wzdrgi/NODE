import sys
from sko.tools import set_run_mode
from sko.DE import DE
from scipy.optimize import minimize
import random
import numpy as np
sys.path.append(r'D:/vspy/.venv/NODE')
import NODE as nd

st_data,sc_data,cell_type,st_pixel = nd.get_test_data(6,50,20,4)

result_data,W_interaction = nd.get_deconvolution(st_data = st_data,
                                    sc_data = sc_data,
                                    cell_type = cell_type,
                                    st_coordinate = st_pixel,
                                    method_optimize = 1,
                                    prossecing_reserve = True,
                                    file_path = 'C:/Users/HP/Desktop/test_file/test_result.txt',
                                    Number_of_iterations = 500)

print(result_data)
print(W_interaction)

