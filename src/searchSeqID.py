import requests
from bs4 import BeautifulSoup
from queue import Queue
from threading import Thread
import numpy as np
import time


### ip池
def get_ip_list():
	"""获取代理ip池
	"""
	url = "https://www.xicidaili.com/wn"
	headers = {
		'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.110 Safari/537.36'
		}
	res = requests.get(url, headers=headers)
	html = res.content.decode()

	soup = BeautifulSoup(html, 'html.parser')
	trs = soup.find_all('tr')[1:]
	ip_list = []
	for tr in trs:
		tds = [td.text for td in tr.find_all('td')]
		ip_addr = tds[1]
		port = tds[2]
		ip_type = tds[5]
		ip_list.append("%s://%s:%s" %(ip_type, ip_addr, port))
	return ip_list


def random_choose_ip(ip_list):
	"""从代理ip池中随机选择ip
	"""
	ip = np.random.choice(ip_list)
	return ip


def searchSeqID(seq, threshold=0.95):
	"""搜索基因序列并返回匹配的 sequence ID 
		params:
			seq: (str) 基因序列
			threshold: (float) 搜索的基因匹配程度,范围[0, 1],
		return:
			(array) 
			[[seq_id1, score1],
			 [seq_id2, score2],
			 ...,
			 [seq_idn, scoren]]
	"""
	sess = requests.Session()
	headers = {
		"User-Agent":"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/71.0.3578.98 Safari/537.36"
	}
	sess.headers.update(headers)
	try:
		ip = random_choose_ip(ip_list)
	except:
		ip_list = get_ip_list()
		ip = random_choose_ip(ip_list)
	proxies = {
		ip.split(':')[0]:ip,
	}
	data = {
		'printParams': 'no',
		'fromDB': 'no',
		'sequence_file': '(binary)',
		'sequence': "%s" %(seq),
		'strain': 'both',
		'source': 'isolates',
		'size': 'gt1200',
		'quality': 'good',
		'taxonomy': 'rdpHome',
		'num': 20,
		'submit': 'Submit',
	}
	url = "http://rdp.cme.msu.edu/seqmatch/SeqmatchControllerServlet/start"
	res = sess.post(url, data=data, proxies=proxies)
	print("已发送数据", res, ip)

	# 必须 Status --> Summary --> Result 的顺序依次访问，否则可能会报错
	# Status
	res = sess.get("http://rdp.cme.msu.edu/seqmatch/seqmatch_status.jsp?qvector=12&depth=0&currentRoot=0&num=20",
					proxies=proxies)
	print("已获取 Status 页面")
	# Summary
	res = sess.get("http://rdp.cme.msu.edu/seqmatch/seqmatch_sum.jsp?qvector=12&depth=0&currentRoot=0&num=20",
					proxies=proxies)
	print("已获取 Summaray 页面")
	# Result
	res = sess.get('http://rdp.cme.msu.edu/seqmatch/seqmatch_result.jsp?qvector=204&depth=0&currentRoot=0&num=20',
				    proxies=proxies)
	html = res.content.decode()
	print("已获取序列比对数据", res)
	
	# 搜索 match_score 大于等于0.95的序列id
	soup = BeautifulSoup(html, 'html.parser')
	detail_html = soup.find('div', {'class':'details'})
	seq_id_list = [item.text for item in detail_html.find_all('a', {'target':'_blank'})]
	match_score = np.array([float(item.text) for item in detail_html.find_all('span', {'style':'background-color: #FFDBB8'})[1:]])
	seq_id_searched = []
	if True not in (match_score >= threshold):
		print('match_score 全部小于 %.2f' %(threshold))
	for i, score in enumerate(match_score):
		if score >= threshold:
			seq_id_searched.append([seq_id_list[i], score])
	return seq_id_searched


def searchBatchSeqID(strain_name, fasta_list, queue):
	"""批次任务
	"""
	batch_data = [strain_name, []]
	for i, seq in enumerate(fasta_list):
		time.sleep(np.random.uniform(0, 3))
		print("%d " %(i+1), end='')
		try:
			seq_id_searched = searchSeqID(seq)
		except:
			seq_id_searched = []
		batch_data[1].extend(seq_id_searched)
	queue.put(batch_data)
	print(batch_data)
	print("成功获取 %s 的所有 seq_id" %(strain_name))


def searchAllSeqID(fasta_need_search):
	"""多线程获取全部菌属的匹配的 seq_id
	"""
	queue = Queue()
	thds = []
	for strain_name, fasta_list in fasta_need_search.items():
		thd = Thread(target=searchBatchSeqID,
					args=(strain_name, fasta_list, queue))
		thd.start()
		thds.append(thd)
		print("开启 %s 的线程" %(strain_name))
	for thd in thds:
		thd.join()
		
	seq_id_searched = []
	for thd in thds:
		seq_id_searched.append(queue.get())
	print('成功获取所有seq_id')
	return seq_id_searched


def main(folder):
	# 读取 fasta_need_search.txt 并获得需要搜索的菌属名与序列信息
	fpath = './dataset/%s/fasta_need_search.txt' %(folder)
	with open(fpath, 'r') as f:
		fasta_need_search = eval(f.read())

	# A test
	# fasta = {"k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;Other": 
	# ["TACGGAGGgggTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCgcgcgTAGGCGGACTAGTCAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGCTGGTCTTGAGTTCGAgagagGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGCCAGTCGTCGGGCAGTATACTGTTCGGTGACacacCTAACGGATTAAGCATTCCGCCTG",
	# "TACGGAGGgggCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACTATTAAGTCAGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTTTGATACTGGTAGTCTAGAGTTCGAgagagGTGAGTGGAACTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGCCAGACGTCGGCAAGCATGCTTGTCGGTGTCACACCTAACGATTAAGCATTCCGGCCTGGGGAGTACGGTCGCAAGATTA"]}
	# seq_id_searched = searchAllSeqID(fasta)
	
	seq_id_searched = searchAllSeqID(fasta_need_search)
	with open('./dataset/%s/seq_id_searched.txt' %(folder), 'w') as f:
		f.write(str(seq_id_searched))
	f.close()


if __name__=="__main__":
	ip_list = get_ip_list()
	main(input("请输入数据集名: "))
