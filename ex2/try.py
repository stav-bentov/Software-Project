import pandas as pd

data1= pd.read_csv("input_1_db_1.txt", header=None) # Put data from file_1 to a data1
data2= pd.read_csv("input_1_db_2.txt", header=None) # Put data from file_2 to a data2
mergeDf=pd.merge(data1, data2, on=0)
mergeDf = mergeDf.sort_values(by=[0])
print(mergeDf)
data_points_array = mergeDf.to_numpy()  # Saves the data points
Index_array = data_points_array[:, 0].astype(int)
print(Index_array)