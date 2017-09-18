import re

def time_of_file(file_name):
    # idx1 = file_name.index('(')
    # idx2 = file_name.index(')')
    # idx3 = file_name.index('-',idx1,idx2)
    # idx4 = file_name.rindex('-',idx1,idx2)
    # h = int(file_name[idx1+1:idx3])
    # m = int(file_name[idx3+1:idx4])
    # s = int(file_name[idx4+1:idx2])
    # return (h,m,s)
    m = re.findall('.*\((\d+)-(\d+)-(\d+).*', file_name)
    return tuple(map(int, m[0]))

def force_of_file(file_name):
    m = re.findall('.*_(\d+(\.\d+|))kn.*', file_name)
    return float(m[0][0])