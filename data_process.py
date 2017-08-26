#-*- coding:utf-8 -*-
import pandas as pd
import sys

cols = ['时间','振次','轴位移','力','应力','变形',
        '侧位移','围压','孔压','反压','体变','小体变']
def data_process(f):
    df = pd.read_table(f, encoding='gb2312')
    df.columns = cols
    df.drop(['侧位移','体变','小体变','应力'],1,inplace=True)
    df.dropna(inplace=True)
    df.to_csv(sys.argv[1][:-4]+'.csv',index=False)

if __name__ == '__main__':
    for file in sys.argv:
        data_process(file)
