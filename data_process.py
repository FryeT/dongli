#-*- coding:utf-8 -*-
import pandas as pd
import numpy as np

cols = ['时间', '振次', '轴位移', '力', '应力', '变形',
        '侧位移', '围压', '孔压', '反压', '体变', '小体变']
drop_cols = ['侧位移', '体变', '小体变', '应力']
A0 = np.pi * 0.3**2 / 4    # area, unit:m^2
H0 = 600                   # height, unit:mm
def data_process(fname):
    """
    把txt格式数据经过简单处理后保存成csv格式
    """
    df = pd.read_table(fname)
    df.columns = cols
    df.drop(drop_cols, 1, inplace=True) #remove cols that not needed
    df.dropna(inplace=True)
    df['sigma_d'] = (df['力']-df['力'][0])/A0
    df['epsilon_d'] = (df['变形']-df['变形'][0])/H0
    df['d1'] = (df['轴位移']-df['轴位移'][0])/H0
    df['pore_pressure'] = df['孔压']-df['孔压'][0]
    df['back_pressure'] = df['反压']-df['反压'][0]
    df.to_csv(fname[:-4]+'.csv', index=False)

if __name__ == '__main__':
    import sys
    for file in sys.argv[1:]:
        data_process(file)
