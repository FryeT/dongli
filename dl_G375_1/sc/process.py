import pandas as pd

df = pd.read_table('datas/dl_G375_', encoding='gb2312')
df.columns=['时间','振次','轴位移','力','应力','变形','侧位移','围压','孔压','反压','体变','小体变']
df.drop(['侧位移','体变','小体变','应力'],1,inplace=True)
df.dropna(inplace=True)

plt.plot(df['变形'],df['力'])