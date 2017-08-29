def time_of_file(file_name):
    idx1 = file_name.index('(')
    idx2 = file_name.index(')')
    idx3 = file_name.index('-',idx1,idx2)
    idx4 = file_name.rindex('-',idx1,idx2)
    h = int(file_name[idx1+1:idx3])
    m = int(file_name[idx3+1:idx4])
    s = int(file_name[idx4+1:idx2])
    return (h,m,s)

