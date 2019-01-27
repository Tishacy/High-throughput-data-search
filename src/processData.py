import pandas as pd
import numpy as np
import os
import re
import json


def needSearch(strain_item):
	"""
		判断菌属是否需要查
		返回值：(need_search, code)
			need_search: True/False 如果需要查则返回True,反之是False
			code: 0/1/2/3 标识码
				0 需要查
				1 菌属信息不全
				2 f 中含有 unidentified
				3 数值全部小于0.005（0.5%）
	"""
	# 解析菌属字段
	strain_name = strain_item.name
	strain_values = strain_item.values

	# 将菌属字符串分割成为k,p,c,o,f,g的列表
	strain_detail = strain_name.split(';')
	
	need_search = True
	code = 0
	# 检查是否有菌属信息
	for detail in strain_detail[:-1]:
		if detail[-1] == '_' or 'Other' in detail:
			need_search = False
			code = 1
			break

	# 检查 f 是否有 ‘unidentified’
	if 'unidentified' in strain_detail[4]:
		need_search = False
		code = 2

	# 检查数值是否小于 0.5%
	if True not in (strain_values > 0.005):
		need_search = False
		code = 3
	
	return need_search, code


def strainNeedSearch(strain_df):
    """
        返回所有需要查找的菌属
    """
    strain_need_search = []
    for i in range(1, len(strain_df)):
        strain_item = strain_df.iloc[i]
        if needSearch(strain_item)[0] == True:
            # 将菌属名存入 strain_name
            # strain_name = ';'.join(strain_item.name.split(';')[:-1])
            strain_name = strain_item.name
            strain_need_search.append(strain_name)
    # print("共需要查 %d 个菌属" %(len(strain_need_search)))
    return strain_need_search


def OTUsNeedSearch(strain_need_search, OTUs_content):
	"""
		根据 strain_need_search 找出对应的 OTUs，
		并存入 OTUs_need_search 字典中
	"""
	OTUs_need_search = dict()

	for strain_name in strain_need_search:
		# 正则表达式找出菌属名所对应的 OTUs_content
		info = r'\n([^\n]+%s[^\n]+)\n' %(';'.join(strain_name.split(';')[:-1]))
		OTUs_content_need_search = re.findall(info, OTUs_content)

		# 从 OTUs_content 中获取 OTUs 列表
		OTUs_list_need_search = [content.split('\t')[0] for content in OTUs_content_need_search]

		# 将 OTUs 列表存入 OTUs_need_search 字典中
		OTUs_need_search[strain_name] = OTUs_list_need_search
	
	return OTUs_need_search


def fastaNeedSearch(OTUs_need_search, fasta_content):
	"""
		根据 OTUs_need_search 找出对应的基因序列，
		并存入 fasta_need_search 字典中
	"""
	fasta_need_search = dict()
	
	for strain_name, OTUs_list in OTUs_need_search.items():
		# 正则表达式找出 OTU 所对应的 fasta_content
		fasta_content_need_search = []
		for OTU in OTUs_list:
			info = r'(%s[^>]+)>' %(OTU)
			fasta_content_need_search.append(re.findall(info, fasta_content)[0])
		
		# 从 fasta_content 中获取基因序列
		fasta_list_need_search = [''.join(content.split('\n')[1:]) for content in fasta_content_need_search]
		
		# 将基因序列存入 fasta_need_search 字典中
		fasta_need_search[strain_name] = fasta_list_need_search
	
	return fasta_need_search


def saveOTUsNeedSearch(OTUs_need_search, folder):
	"""
		将 OTUs_need_search 保存至json文件
	"""
	with open('./dataset/%s/OTUs_need_search.txt' %(folder),'w') as outfile:
		# json.dump(OTUs_need_search, outfile)
		outfile.write(str(OTUs_need_search))
	outfile.close()


def saveFastaNeedSearch(fasta_need_search, folder):
	"""
		将 fasta_need_search 保存至json文件
	"""
	with open('./dataset/%s/fasta_need_search.txt' %(folder),'w') as outfile:
		# json.dump(fasta_need_search, outfile)       
		outfile.write(str(fasta_need_search))
	outfile.close()


def main(folder):
	# 打开数据文件
	strain_df = pd.read_excel('./dataset/%s/strain.xlsx' %(folder))
	with open('./dataset/%s/OTUs.tax_assignments.txt' %(folder), 'r') as f:
		OTUs_content = f.read()
	with open('./dataset/%s/OTUs.fasta' %(folder), 'r') as f:
		fasta_content = f.read()

	# 执行操作
	strain_need_search = strainNeedSearch(strain_df)
	print('已找到所有待搜索菌属，共 %d 个菌属' %(len(strain_need_search)))
	OTUs_need_search = OTUsNeedSearch(strain_need_search, OTUs_content)
	print('已找到共 %d 个待搜索菌属的OTU数' %(len(OTUs_need_search)))
	fasta_need_search = fastaNeedSearch(OTUs_need_search, fasta_content)
	print('已找到共 %d 个待搜索菌属所对应的基因序列' %(len(fasta_need_search)))
	saveOTUsNeedSearch(OTUs_need_search, folder)
	saveFastaNeedSearch(fasta_need_search, folder)
	print('已保存数据文件')


if __name__=="__main__":
	main(input('请输入数据集名: '))
