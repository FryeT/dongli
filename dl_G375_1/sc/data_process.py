import pandas as pd
import sys

df = pd.read_table(sys.argv[1], encoding='gb2312')
df.columns=['时间','振次','轴位移','力','应力','变形','侧位移','围压','孔压','反压','体变','小体变']
df.drop(['侧位移','体变','小体变','应力'],1,inplace=True)
df.dropna(inplace=True)
df.to_csv(sys.argv[1][:-4]+'.csv',index=False)
