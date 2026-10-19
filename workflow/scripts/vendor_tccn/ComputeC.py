import numpy as np
import itertools 

def is_subarray(arr, subarr):
    """
    Check if subarr is a subarray of arr.

    Parameters:
    arr (list): The main array.
    subarr (list): The subarray to check.

    Returns:
    bool: True if subarr is a subarray of arr, False otherwise.
    """
    # Get the lengths of the main array and the subarray
    n, m = len(arr), len(subarr)
    
    # If the subarray length is greater than the main array length, it cannot be a subarray
    if m > n:
        return False

    # Traverse the main array to check for the subarray
    for i in range(n - m + 1):
        # Slice the main array to get a part of the same length as the subarray
        current_subarr = arr[i:i + m]
        
        # Check if the sliced part is equal to the subarray
        if current_subarr == subarr:
            return True

    # If no matching subarray is found, return False
    return False

def generate_combinations_matrix(num_letters):
    """
    Generate a matrix of all possible combinations of letters.

    Parameters:
    num_letters (int): The number of letters to generate combinations for.

    Returns:
    list: A matrix where each row represents a combination of letters.
    """
    letters = [chr(ord('A') + i) for i in range(num_letters)]  # Automatically generate a list of letters
    n = len(letters)
    combinations = []

    # Generate all possible combinations
    for r1 in range(1, n + 1):
        for combo in itertools.combinations(letters, r1):
            row = [1 if letter in combo else 0 for letter in letters]
            combinations.append(row)

    # Sort by the number of 1s in ascending order
    combinations.sort(key=lambda row: sum(row))

    return combinations

def find_rows_with_one_more(matrix, row):
    """
    Find rows in the matrix that have exactly one more '1' than the given row.

    Parameters:
    matrix (list of list of int): The input 0-1 matrix.
    row (int): The index of the reference row.

    Returns:
    list: A list of row indices that have exactly one more '1' than the reference row.
    """
    reference_row = matrix[row]
    linenumber = []
    columns_of_ones = []

    # Find the indices of '1's in each row
    for row2 in matrix:
        ones_indices = [i for i, value in enumerate(row2) if value == 1]
        columns_of_ones.append(ones_indices)

    # Find rows that have exactly one more '1' and the same positions for the original '1's
    for row3, j in enumerate(matrix):
        if np.sum(j) == np.sum(reference_row) + 1 and set(columns_of_ones[row]).issubset(columns_of_ones[row3]):
            linenumber.append(row3)
    return linenumber

def find_rows_by_indices(matrix, indices):
    """
    Find rows in the matrix by their indices.

    Parameters:
    matrix (list of list of int): The input matrix.
    indices (list of int): The list of row indices to select.

    Returns:
    list: A list of selected rows from the matrix.
    """
    selected_rows = [matrix[i] for i in indices]
    return selected_rows

def find_int(num):
    """
    Custom integer rounding function.

    Parameters:
    num (float): The input number.

    Returns:
    int: The rounded integer.
    """
    if int(num) == num:
        return int(num) - 1
    else:
        return int(num)
    
import numpy as np

def computeC(arrays_gene, readlengh, gene_dict):
    """
    Compute the C matrix based on the given arrays, read length, and gene dictionary.

    Parameters:
    arrays_gene (list of list of str): The input arrays of gene identifiers.
    readlengh (int): The read length.
    gene_dict (dict): A dictionary mapping gene identifiers to their lengths.

    Returns:
    np.ndarray: The computed C matrix.
    """
    # Create a copy of arrays
    arrays = [array[:] for array in arrays_gene]
    
    # Remove arrays with length less than readlengh
    length = []
    for array in arrays:
        length_array = sum([gene_dict[key] for key in array])
        length.append(length_array)
    
    # Check if the maximum length is greater than readlengh
    if max(length) < readlengh:
        return 'The read length is too long, please try a smaller one.'

    # Filter arrays based on length
    bools = [x > readlengh for x in length]
    arrays = [x for x, keep in zip(arrays, bools) if keep]
    
    # Generate the r matrix
    r = generate_combinations_matrix(len(arrays))
    r = [i for i in r if i.count(1) > 1]

    # Find the common parts for each combination
    result_array = []
    for i in r:
        k = 0
        for j in range(len(arrays)):
            if i[j] == 1 and k == 0:
                common_result = set(arrays[j])
                k = 1
            elif i[j] == 1 and k == 1:
                common_result &= set(arrays[j])
        result_array.append(list(common_result))
    result_array = [sorted(i) for i in result_array]

    # Match the common elements with the r matrix
    result_array_part = []
    for index, row in enumerate(result_array, 1):
        result = []
        location = []
        location_array = []
        if len(row) == 0:
            result = []
        else:
            current_group = [row[0]]
            k = 0
            for num in row[1:]:
                # Check which transcripts are used in the r matrix
                for i in range(len(arrays)):
                    if r[index-1][i] == 1:
                        location.append(arrays[i].index(num))
                        location_array.append(i)
                        
                # Determine if the previous element is connected in the transcripts
                for j in range(len(location_array)):
                    if arrays[location_array[j]][location[j]-1] == row[row.index(num) - 1]:
                        k += 1
                if k == len(location_array):
                    current_group.append(num)
                else:
                    result.append(current_group)
                    current_group = [num]
                location_array = []
                location = []
                k = 0
            result.append(current_group)
        result_array_part.append(result)
        
    # Remove elements with length less than readlengh
    realans = []
    for i in range(len(result_array_part)):
        ans = [sublist for sublist in result_array_part[i] if sum([gene_dict[key] for key in sublist]) >= readlengh]
        realans.append(ans)
    
    # Calculate the length of the common parts
    totalsum = []
    for row, overlapping in enumerate(realans):
        # Find rows with one more '1', denoted as overpart
        overpart = find_rows_by_indices(realans, find_rows_with_one_more(r, row))
        
        # Calculate the last row separately
        if overpart == []:
            npublic = 0
            for i in realans[len(realans)-1]:
                npublic = (sum([gene_dict[key] for key in i]) - readlengh + 1) + npublic
            totalsum.append(npublic)
        else:
            # Calculate the overlapping parts
            sum_result2 = 0
            for i in overlapping:
                iii = i

                # Find the parts with one more '1' and remove duplicates
                includeii = []
                for rowii in overpart:
                    for elementii in rowii:
                        if is_subarray(iii, elementii):
                            includeii.append(elementii)
                            
                # Remove duplicates
                includeii = [list(x) for x in set(tuple(x) for x in includeii)]
                
                # Remove fully contained elements
                includeiinew = []
                for element2 in includeii:
                    suminclude = 0
                    for j in includeii:
                        if is_subarray(j, element2):
                            suminclude += 1
                    if suminclude == 1:
                        includeiinew.append(element2)
                includeiinew = sorted(includeiinew, key=lambda x: x[0])

                # Process includeiinew
                if len(includeiinew) == 0:
                    sum_result2 += sum([gene_dict[key] for key in i]) - readlengh + 1
                else:
                    # Process i in overlapping
                    sum_i = sum([gene_dict[key] for key in i])
                    
                    # Generate a corresponding 0 array
                    list_i = [0 for _ in range(sum_i - readlengh + 1)]

                    # Process elements in includeiinew
                    for element_i in includeiinew:
                        # Find the first element and its position in i
                        first_i = element_i[0]
                        location_i = i.index(first_i)
                        
                        # Calculate the length of element_i
                        sum_ele_i = sum([gene_dict[key] for key in element_i])
                        
                        # Update list_i based on the position and length of element_i
                        sum_first = sum([gene_dict[key] for key in i[:location_i]])  # Calculate the length before element_i
                        list_i[sum_first:sum_ele_i + sum_first - readlengh + 1] = [1] * (sum_ele_i - readlengh + 1)

                    # Count the number of 0s in list_i
                    count_of_zeros = list_i.count(0)
                    sum_result2 += count_of_zeros
            totalsum.append(sum_result2)
    
    # Calculate the C matrix
    n = len(arrays)
    rows = 2**n - 1
    cols = n
    Cmatrix = np.zeros((rows, cols))
    for i in range(n, 2**n-1):
        for j in range(n):
            if r[i-n][j] == 1:
                Cmatrix[i][j] = totalsum[i-n] / (sum([gene_dict[key] for key in arrays[j]]) - readlengh + 1)
                
    # Calculate the sum of each column
    column_sums = np.sum(Cmatrix, axis=0)
    for k in range(n):
        Cmatrix[k][k] = 1 - column_sums[k]
        
    return Cmatrix