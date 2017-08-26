#-*- coding:utf-8 -*-

def modify(f, d):
    with open(f,'r+',encoding='gb2312') as f:
        s = f.read().replace(d,'\t')
        s = s.lstrip()
        f.seek(0)
        f.write(s)

if __name__ == '__main__':
    from sys import argv
    import os
    for f in argv[1:]:
        if os.path.isfile(f):
            modify(f,',')


