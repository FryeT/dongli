import glob

def modify(files, d):
    for file in files:
        with open(file,'r+') as f:
            s = f.read().replace(d,'\t')
            s = s.lstrip()
            f.seek(0)
            f.write(s)

files = glob.glob('*.txt')
modify(files, ',')

