# NODE
A Python program for deconvolution spatial transcriptomics data and inference of spatial communication
# Background
NODE is a Python-based package for spatial transcriptomic data. NODE is based on an optimization search model and non-negative least squares problem to use scRNA-Seq data for deconvoluting spatial transcriptomics data and inferring spatial communications. In deconvolution, Node can infer cell number and cell type from spatial transcriptomics data by referring to single-cell data. In inference of spatial communication, NODE can infer the information flow in space.
# Usage
NODE require users to provide four types of data for deconvolution and inferring the spatial communication.
The specific form of the data is shown below:
## Sc_data (Single-cell data) 
Sc_data needs to carry the gene name and cell name as a square, (0, 0) can specify any name, such as 'gene' and 'name'. Sc_data rows represent genes, columns represent cells, and each column represents the expression of a gene. The genes' name of sc_data must correspond to the genes' name of st_data.

sc_data is shown below:

    name      cell1      cell2 

    gene1     5.0        3.0

    gene2     2.0        0.0

## St_data (Spatial transcriptomics data) 
St_data needs to carry the genes' name and spots' name as a square,(0,0) can specify any name, such as 'spot' and 'name'. st_data rows represent genes, columns represent spots, and each column represents the expression of a gene in a spot. The genes' name of st_data must correspond to the genes's name of sc_data.

st_data is shown below:

    spot     spot1     spot2

    gene1     4.0      8.0

    gene2     0.0      2.0

## cell_type (single-cell categorization data)
cell_type needs to carry the cell name and its cell type. celltype's first line is title information including name and cell type name. Cell_type rows represent each cell. The name column of cell_type must correspond to the first row of sc_data, and the cell names in it must be the same.

cell_type is shown below:
    
    name          celltype

    cell1         type1

    cell2         type2

## st_coordinate (coordinate data of spatial transcriptomic data)
st_coordinate needs to carry the spots' name and theirs coordinates. St_coordinate 's first line is title information including 'spot', 'x', and 'y'. St_coordinate rows represent spot. St_coordinate's spot column must correspond to the first row of st_data, where the names of the spots must be the same.

st_coordinate is shown below:

    spot    x       y

    spot1   8       20

    spot2   10      12
## get_deconvolution
Once the user has prepared the above four sets of data, the deconvolution can be performed with the get_deconvolution function. 
The method_optimize parameter has two options, 1 and 2, with method 1 being faster and method 2 having tighter control over conditions.

    from sko.tools import set_run_mode
    from sko.DE import DE
    from scipy.optimize import minimize
    import random
    import numpy 
    import sys
    sys.path.append(r'/path/')
    #The path is the location of the user's NODE program.
    import NODE as nd
    
    st_data,sc_data,cell_type,st_coordinate = nd.get_test_data(spot_number = 6, #spot's number
                                                            gene_number = 50,#gen's number
                                                            cell_number = 20,#cell's number
                                                            type_number = 4  #cell type's number
                                                            )
    # This is only used to generate test data, the content of the data is meaningless,
    
    result_data,spatial_communication_result = nd.get_deconvolution(    
                                                 st_data = st_data,
                                                 sc_data = sc_data,
                                                 cell_type = cell_type,
                                                 st_coordinate = st_coordinate,
                                                 neighbors_number = 4,
                                                 Number_of_iterations = 500,
                                                 method_optimize = 1)
                                                 #result_data is deconvolution result
                                                 #spatial_communication_result is the spatial communication result
## deconvolution result
    use scipy to deconvolution
    Deconvolution Start
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 188.74543074614562
                Iterations: 9
                Function evaluations: 49
                Gradient evaluations: 9
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 176.72391283154244
                Iterations: 6
                Function evaluations: 56
                Gradient evaluations: 6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 176.02826151754223
                Iterations: 5
                Function evaluations: 28
                Gradient evaluations: 5
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 211.6447486426664
                Iterations: 10
                Function evaluations: 53
                Gradient evaluations: 10
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 202.85144307047065
                Iterations: 7
                Function evaluations: 61
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 202.18242184095993
                Iterations: 6
                Function evaluations: 34
                Gradient evaluations: 6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 201.94078423039363
                Iterations: 8
                Function evaluations: 45
                Gradient evaluations: 8
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 196.01292063238918
                Iterations: 6
                Function evaluations: 55
                Gradient evaluations: 6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 195.34932818164293
                Iterations: 6
                Function evaluations: 34
                Gradient evaluations: 6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 196.67580968073435
                Iterations: 8
                Function evaluations: 45
                Gradient evaluations: 8
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 188.41134772726363
                Iterations: 7
                Function evaluations: 61
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 188.06598054373694
                Iterations: 8
                Function evaluations: 43
                Gradient evaluations: 8
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 171.98296870090658
                Iterations: 7
                Function evaluations: 42
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 170.01704496370945
                Iterations: 7
                Function evaluations: 60
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 169.84450325348197
                Iterations: 7
                Function evaluations: 39
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 175.05910614978137
                Iterations: 8
                Function evaluations: 44
                Gradient evaluations: 8
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 167.34917126671047
                Iterations: 7
                Function evaluations: 61
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 168.68511851721274
                Iterations: 7
                Function evaluations: 39
                Gradient evaluations: 7
    Deconvolution Finish!
    
## get_deconvolution with save path
The method_optimize parameter has two options, 1 and 2, with method 1 being faster and method 2 having tighter control over conditions.

    from sko.tools import set_run_mode
    from sko.DE import DE
    from scipy.optimize import minimize
    import random
    import numpy 
    import sys
    sys.path.append(r'/path/')
    #The path is the location of the user's NODE program.
    import NODE as nd
    
    st_data,sc_data,cell_type,st_coordinate = nd.get_test_data(spot_number = 6, #spot's number
                                                            gene_number = 50,#gen's number
                                                            cell_number = 20,#cell's number
                                                            type_number = 4  #cell type's number
                                                            )
    # This is only used to generate test data, the content of the data is meaningless,
    
    result_data,spatial_communication_result = nd.get_deconvolution(    
                                                 st_data = st_data,
                                                 sc_data = sc_data,
                                                 cell_type = cell_type,
                                                 st_coordinate = st_coordinate,
                                                 neighbors_number = 4,
                                                 Number_of_iterations = 500,
                                                 method_optimize = 1,
                                                 prossecing_reserve = True,
                                                 file_path = "/save path/")
                                                 #result_data is deconvolution result
                                                 #spatial_communication_result is the spatial communication result
                                                 #The "save path" is the path where the user wants to store the result of the deconvolution
                                                 
## deconvolution result
    use scipy to deconvolution
    Deconvolution Start!
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 196.72161593468078
                Iterations: 9
                Function evaluations: 49
                Gradient evaluations: 9
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 188.6195490890819
                Iterations: 7
                Function evaluations: 60
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 188.1374912940327
                Iterations: 6
                Function evaluations: 34
                Gradient evaluations: 6
    Saving during the process1/6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 215.5801982094186
                Iterations: 8
                Function evaluations: 45
                Gradient evaluations: 8
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 204.2401902722098
                Iterations: 7
                Function evaluations: 60
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 203.7423141369403
                Iterations: 6
                Function evaluations: 34
                Gradient evaluations: 6
    Saving during the process2/6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 205.3893501761546
                Iterations: 8
                Function evaluations: 45
                Gradient evaluations: 8
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 197.02788900931657
                Iterations: 7
                Function evaluations: 60
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 196.0322029957753
                Iterations: 6
                Function evaluations: 34
                Gradient evaluations: 6
    Saving during the process3/6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 219.7033689902071
                Iterations: 9
                Function evaluations: 49
                Gradient evaluations: 9
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 214.41871359503847
                Iterations: 7
                Function evaluations: 59
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 213.10491917916045
                Iterations: 7
                Function evaluations: 38
                Gradient evaluations: 7
    Saving during the process4/6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 185.62308511046183
                Iterations: 9
                Function evaluations: 50
                Gradient evaluations: 9
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 181.30355286463137
                Iterations: 7
                Function evaluations: 60
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 180.39147280274446
                Iterations: 6
                Function evaluations: 34
                Gradient evaluations: 6
    Saving during the process5/6
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 174.7447600239144
                Iterations: 10
                Function evaluations: 55
                Gradient evaluations: 10
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 167.0204664221119
                Iterations: 7
                Function evaluations: 60
                Gradient evaluations: 7
    Optimization terminated successfully    (Exit mode 0)
                Current function value: 166.12674797324232
                Iterations: 6
                Function evaluations: 35
                Gradient evaluations: 6
    Saving during the process6/6
    Deconvolution Finish!
