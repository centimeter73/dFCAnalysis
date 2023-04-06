'''
参考施雪蓉组会分享的文章，进行时变性分析
Abnormal dynamic functional connectivity after sleep deprivation from temporal variability perspective
'''

import os
import numpy as np
from statsmodels.stats import multitest
from statsmodels.stats.multitest import multipletests
from scipy import stats
from scipy.stats import ttest_ind
import statsmodels


def statsAnalysis(regional_variability_values_all,node_number,DM_number,HC_number):
	p_values = []
	for t in range(0,node_number):
		data_all = []
		for row in range(0,len(regional_variability_values_all)):
			data_all.append(regional_variability_values_all[row][t])
		data1 = data_all[0:DM_number]
		data2 = data_all[DM_number:]

		# 双样本t检验
		p_value = stats.ttest_ind(np.array(data1),np.array(data2))[1]
		p_values.append(p_value)


	# 多重比较校正(bh校正)
	p_values_corrected = statsmodels.stats.multitest.multipletests(p_values,alpha = 0.05,method = 'fdr_bh')[1]

	flag = False
	for p_value_corrected in p_values_corrected:
		if p_value_corrected < 0.05:
			flag = True

	if flag:
		print("存在差异节点")

	else:
		print("无显著性差异结果")


def mainAnalysis(root_path,node_number):
	regional_variability_values_all = []
	for folder in os.listdir(root_path):
		subject_path = root_path + "/" + folder
		data_all = {} # 单被试的所有滑动窗功能连接矩阵
		count = 0
		for file in os.listdir(subject_path):
			count = count + 1
			data = np.loadtxt(subject_path + "/" + file)
			data_all[str(count)] = data

		regional_variability_values = []
		regional_variability_matrix = []
		for t in range(0,node_number):
			for value in data_all.values():
				regional_variability_matrix.append(value[t])
			corr_matrix = np.corrcoef(regional_variability_matrix,rowvar = 1)
			corr_values = []
			for p in range(0,len(corr_matrix)):
				for q in range(0,len(corr_matrix[p])):
					if p < q:
						corr_values.append(corr_matrix[p][q])
			corr_avg = np.mean(corr_values)
			regional_variability_values.append(1 - corr_avg)
		regional_variability_values_all.append(regional_variability_values)
	return regional_variability_values_all


# 主函数
if __name__ == '__main__':
	root_path = "I:/T2DM/rsfMRI_FunctionalNetworks/rsfMRI_DFC_latest/246_50TR_120_nonoverlapping"
	node_number = 246
	regional_variability_values_all = mainAnalysis(root_path,node_number)
	# 统计分析
	DM_number = 68
	HC_number = 52
	statsAnalysis(regional_variability_values_all,node_number,DM_number,HC_number)
