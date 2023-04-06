'''
1.基于滑动窗口法构建单被试的系列动态脑功能连接矩阵(z-transformed)
2.基于方差概念计算功能连接网络的时间变异性
3.组间统计分析
考虑问题：稀疏化怎么做(暂时不考虑)
负连接如何处理？ 一般是直接置0
Shared and distinct patterns of dynamic functional connectivity variability of thalamo-cortical circuit in
bipolar depression and major depressive disorder
'''
import os
import numpy as np
from statsmodels.stats import multitest
from statsmodels.stats.multitest import multipletests
from scipy import stats
from scipy.stats import ttest_ind
import statsmodels
	
def mainAnalysis(root_path,node_number):
	for folder in os.listdir(root_path):
		# 窗口矩阵整合
		subject_path = root_path + "/" + folder
		subject_info = {}
		count = 0
		for file in os.listdir(subject_path):
			count = count + 1
			dynamic_functional_matrix = []
			with open(subject_path + "/" + file) as f:
				for line in f.readlines():
					row = line.strip().split("\t")
					new_row = []
					for e in row:
						new_row.append(eval(e))
					dynamic_functional_matrix.append(new_row)
			subject_info[str(count)] = dynamic_functional_matrix

		# 方差时变分析
		variance_matrix = []
		for p in range(0,node_number):
			variance_row = []
			for q in range(0,node_number):
				node_values = []
				for matrix in subject_info.values():
					node_values.append(matrix[p][q])
				node_variance = np.var(node_values)
				variance_row.append(node_variance)
			variance_matrix.append(variance_row)

		# 结果文件输出
		np.savetxt(result_path + "/" + folder + ".txt",variance_matrix,delimiter = "\t",fmt = "%f")


def connAnalysis(result_path,node_number,DM_number,HC_number):
	conns_all = []
	for file in os.listdir(root_path):
		data = np.loadtxt(root_path + "/" + file)
		
		# 每条连接的差异性检验
		conns = []
		for i in range(0,node_number):
			for j in range(0,node_number):
				if i < j:
					conns.append(data[i][j])
		conns_all.append(conns)
	
	# 双样本t检验
	p_values = []
	for t in range(0,len(conns_all[0])):
		p_value = stats.ttest_ind(np.array(conns_all)[0:DM_number,t],np.array(conns_all)[DM_number:DM_number + HC_number,t])[1]	
		p_values.append(p_value)
	

	# 多重比较校正(bh校正)
	p_values_corrected = statsmodels.stats.multitest.multipletests(p_values,alpha = 0.05,method = 'fdr_bh')[1]

	flag = False
	for p_value_corrected in p_values_corrected:
		if p_value_corrected < 0.05:
			flag = True

	if flag:
		print("存在差异连接")
	else:
		print("无显著性差异结果")



if __name__ == '__main__':
	root_path = "I:/T2DM/rsfMRI_FunctionalNetworks/rsfMRI_DFC_latest/246_50TR_1TR_120_modified"
	result_path = "I:/T2DM/rsfMRI_FunctionalNetworks/rsfMRI_DFC_latest/246_50TR_1TR_120_variance_matrix"
	node_number = 246
	DM_number = 68
	HC_number = 52
	mainAnalysis(root_path,node_number)
	connAnalysis(result_path,node_number,DM_number,HC_number)
