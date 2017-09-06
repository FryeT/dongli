#-*- coding:utf-8 -*-
import pandas as pd
import sys

cols = ['时间','振次','轴位移','力','应力','变形',
        '侧位移','围压','孔压','反压','体变','小体变']
drop_cols = ['侧位移','体变','小体变','应力']
def data_process(fname):
    df = pd.read_table(fname)
    df.columns = cols
    df.drop(drop_cols,1,inplace=True)
    df.dropna(inplace=True)
    df.to_csv(fname[:-4]+'.csv',index=False)

if __name__ == '__main__':
    for file in sys.argv[1:]:
        data_process(file)
