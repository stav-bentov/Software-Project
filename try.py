import pandas as pd
import numpy
data1= pd.read_csv('input_1_db_1.txt', header=None)
data2= pd.read_csv('input_1_db_2.txt', header=None)
mergeDf=pd.merge(data1, data2, on=0)
data_points_array=mergeDf.to_numpy()
print(len(data_points_array))
#print(data1)
#print(data2)
#print(mergeData)