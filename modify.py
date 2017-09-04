#-*- coding:utf-8 -*-

def modify(fname, d):
    try:
        f = open(fname, encoding='utf-8')
        s = f.read()
    except:
        f = open(fname, encoding='gb2312')
        s = f.read().replace(d,'\t')
        f.close()
        s = s.lstrip()
        s = s.encode('utf-8')
        f = open(fname,'wb')
        f.write(s)
        f.close()
        print('modified')
    else:
        f.close()
        print('have been modified')

if __name__ == '__main__':
    from sys import argv
    import os
    for fname in argv[1:]:
        print(fname+':')
        if os.path.isfile(fname):
            modify(fname,',')


