'''
参考杨正武师兄文中的时间变异性定义，进行动态脑功能连接网络的时间变异性分析
后续统计分析使用NBS识别异常的连接子网络

'''
import os
import numpy as np


def mainAnalysis(root_path,result_path,node_number):
	for folder in os.listdir(root_path):
		subject_path = root_path + "/" + folder
		data_all = {}
		count = 0
		for file in os.listdir(subject_path):
			count = count + 1
			data = np.loadtxt(subject_path + "/" + file)
			data_all[str(count)] = data

		data_variance_matrix = [] # 时间变异性矩阵
		for i in range(0,node_number):
			data_variance_row = []
			for j in range(0,node_number):
				values = []
				for value in data_all.values():
					values.append(value[i][j])
				values_std = np.std(np.array(values))
				data_variance_row.append(values_std)
			data_variance_matrix.append(data_variance_row)

		os.mkdir(result_path + "/" + folder)
		np.savetxt(result_path + "/" + folder + "/" + folder + ".txt",data_variance_matrix,delimiter = "\t",fmt = "%f")




if __name__ == '__main__':
	root_path = "I:/T2DM/rsfMRI_FunctionalNetworks/rsfMRI_DFC_latest/246_50TR_1TR_120_modified" # z-transformed
	result_path = "I:/T2DM/rsfMRI_FunctionalNetworks/rsfMRI_DFC_latest/246_50TR_1TR_120_variance_yangzhengwu"
	node_number = 246
	mainAnalysis(root_path,result_path,node_number)
