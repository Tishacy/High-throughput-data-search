import pandas as pd
import requests
import re


def get_genus(seq_id):
	"""根据 sequence ID 搜索 genus
	"""
	seq_id_url = "http://rdp.cme.msu.edu/seqmatch/seqmatch_seqrecorddetail.jsp?seqid=%s" %(seq_id)
	res = requests.get(seq_id_url)
	content = res.content.decode()
	organism_info = r"ORGANISM([^*]+).\nREFERENCE"
	genus = re.findall(organism_info, content)[0].split('\n')[2].split('; ')[-1][:-1]
	return genus


def extract_seq_id(folder):
	"""提取指定文件夹中的 seq_id_searched.txt 的数据
	"""
	# 读取数据
	with open('./dataset/%s/seq_id_searched.txt' %(folder)) as f:
		seq_id_searched = eval(f.read())
		f.close()
	# 提取数据
	info_df = {
		"strain_name": [],
		"genus_name": [],
		"seq_id": [],
		"match_score": [],
	}
	for i, info in enumerate(seq_id_searched):
		info_df["strain_name"].append(info[0])
		if info[1] == []:
			info_df['genus_name'].append('Other')
			info_df["seq_id"].append('Other')
			info_df['match_score'].append('Other')
		else:
			best_match = pd.DataFrame(info[1]).sort_values(by=1, ascending=False).values[0]
			info_df["genus_name"].append(get_genus(best_match[0]))
			info_df["seq_id"].append(best_match[0])
			info_df['match_score'].append(best_match[1])
		print("%3d %s" %(i, info[0]))
	return pd.DataFrame(info_df)


def main(folder):
	df = extract_seq_id(folder)
	df.to_excel('./dataset/%s/best_match_strain.xlsx' %(folder))
	print("Finished.")


if __name__=="__main__":
	main(input("请输入数据集名: "))