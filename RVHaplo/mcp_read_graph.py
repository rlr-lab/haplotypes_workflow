import sys
import concurrent.futures
from tqdm import tqdm
import pysam
import numpy as np
import pickle
import multiprocessing as mp
mp.set_start_method("fork")

### load data
file_in=sys.argv[1]
file_snv = sys.argv[2]
cond_pro = float(sys.argv[3])
smallest_snv = int(sys.argv[4])
num_read_1 = int(sys.argv[5])
num_read_2 = int(sys.argv[6])
gap = int(sys.argv[7])
weight_read = float(sys.argv[8])
ovlap_read = int(sys.argv[9])
file_prefix = sys.argv[10]
fre_most_base = float(sys.argv[11])
thread = int(sys.argv[12])
only_snv = sys.argv[13]
sub_graph = int(sys.argv[14])

bamfile = pysam.AlignmentFile(file_in)
ref_name = bamfile.references[0]
f_snv = open(file_snv,'r')
snv_sites = f_snv.readline().split()
snv_sites = np.array([int(i)-1 for i in snv_sites])
f_snv.close()
ref_length = bamfile.get_reference_length(ref_name)

### data initialization
def error_rate(qual):
    er = 10 ** (-qual / 10)
    return er

###  generate the sequence matrix
def call_seq_mat(para):
	aligned_pairs = para[0]
	aligned_seq = para[1]
	seq_name = para[2]
	map_qual = para[3]
	ref_length = para[4]
	snv_sites = para[5]
	align_index = np.array(aligned_pairs)
	align_seq = np.array(list(aligned_seq))
	read_seq = np.array(['-']*ref_length)
	read_seq[list(align_index[:,1])]= align_seq[list(align_index[:,0])]
	map_acc = 1-error_rate(map_qual)
	flag = 0
	if not (read_seq[snv_sites] == '-').all():
		flag = 1
	return seq_name,map_acc,read_seq[snv_sites],flag

print("Generate sequence matrix at candidate SNV sites")

query_list = [(r.get_aligned_pairs(matches_only=True),r.query_sequence,r.query_name,r.mapping_quality,ref_length,snv_sites) for r in bamfile.fetch()]
with concurrent.futures.ProcessPoolExecutor(thread) as executor:
    res = list(tqdm(executor.map(call_seq_mat, query_list), total=len(query_list)))

Map_acc = {}
Reads_all = []
seq_mat = []

for r in res:
	if r[-1]:
		Reads_all.append(r[0])
		seq_mat.append(r[2])
		Map_acc[r[0]] = r[1]

seq_mat = np.array(seq_mat)

#### count the nucleotide
num_snv = len(snv_sites)
num_reads = len(Reads_all)

def count_acgt(para):
	chara = ['A','C','G','T']
	temp = []
	for j in chara:
		temp.append(np.count_nonzero(para==j))
	return temp

print("Count base number at each candidate SNV site")
query_list = [seq_mat[:,i] for i in range(len(snv_sites))]
with concurrent.futures.ProcessPoolExecutor(thread) as executor:
    res = list(tqdm(executor.map(count_acgt, query_list), total=len(query_list)))

Coun = np.array(res)
del res



####  calculate nucleotide frequency and rank
def rank_nucl(Coun):
	res = []
	for i in Coun:
		res.append(list(np.sort(i)/sum(i)))
	return res

nucl_fre = np.array(rank_nucl(Coun))

###   find out the candicate site to be verified >=fre_most_base  
to_be_verified_sites = np.array([i for i in range(num_snv)])
to_be_verified_sites = to_be_verified_sites[nucl_fre[:,-1]>=fre_most_base]

###  calculate the 1-order  probability for the SNV at each candicate site
subsitution_rate = {i:nucl_fre[i,2]/sum(nucl_fre[i,2:]) for i in to_be_verified_sites}

###  reserve the top-2 bases  
chara = ['A','C','G','T']
snv_base = []
major_base = []
for i in range(len(Coun)):
    temp = np.argsort(Coun[i])
    snv_base.append(chara[temp[-2]])
    major_base.append(chara[temp[-1]])

snv_base = np.array(snv_base)
major_base = np.array(major_base)

seq_snv_flag = []
seq_major_flag = []
for i in range(len(seq_mat)):
    temp = seq_mat[i]==snv_base
    seq_snv_flag.append(temp)
    temp = seq_mat[i]==major_base
    seq_major_flag.append(temp)

seq_snv_flag = np.array(seq_snv_flag)
seq_major_flag = np.array(seq_major_flag)

###   calculate the conditional probability
def cond_pro_1_order(c_site,g_site,ssf=seq_snv_flag,smf=seq_major_flag,v_min = num_read_1):
    #c_site: current site
    #g_site: given site
    #ssf: seq_snv_flag
    #smf: seq_major_flag
    current_snv_flag = ssf[:,c_site]
    current_major_flag = smf[:,c_site]
    given_flag = ssf[:,g_site]
    num_snv = np.count_nonzero(current_snv_flag[given_flag])
    num_major = np.count_nonzero(current_major_flag[given_flag])
    denominator  = num_snv + num_major
    numerator = num_snv
    if denominator > v_min:
        pro = numerator/(denominator+1e-16)
    else:
        pro = 0
    return pro

def cond_pro_extend(c_site,g_site,cg_flag,ssf,smf):
    #c_site: current site
    #g_site: given site
    #ssf: seq_snv_flag
    #smf: seq_major_flag
    #current_given_flag
    current_snv_flag = ssf[:,c_site]
    current_major_flag = smf[:,c_site]
    given_flag = np.logical_and(ssf[:,g_site],cg_flag)
    num_snv = np.count_nonzero(current_snv_flag[given_flag])
    num_major = np.count_nonzero(current_major_flag[given_flag])
    denominator = num_snv + num_major
    numerator = num_snv
    pro = numerator/(denominator+1e-16)
    return pro,given_flag,denominator

def verify_snv_site(site_id,to_be_verified_sites=to_be_verified_sites,seq_snv_flag=seq_snv_flag,seq_major_flag=seq_major_flag,v_min_1=num_read_1,v_min_2=num_read_2,hd=gap):
    temp_flag = to_be_verified_sites!=site_id
    given_sites = to_be_verified_sites[temp_flag]
    con_pro = []
    for j in given_sites:
        con_pro.append(cond_pro_1_order(site_id,j))
    con_pro = np.array(con_pro)
    index_sort = np.argsort(-con_pro)
    current_pro = subsitution_rate[site_id]
    current_given_flag = np.array([True]*num_reads)
    for k in index_sort:
        given_site = given_sites[k]
        if abs(snv_sites[given_site]-snv_sites[site_id]) < hd:
            continue
        new_pro,new_given_flag,new_v = cond_pro_extend(site_id,given_site,current_given_flag,seq_snv_flag,seq_major_flag)
        if new_pro>current_pro and new_v>=v_min_2:
            current_pro = new_pro
            current_given_flag = new_given_flag
            if current_pro>=cond_pro:
                break
        else:
            break
    return current_pro

print("maximum conditional probability")
query_list = [i for i in to_be_verified_sites]

with concurrent.futures.ProcessPoolExecutor(thread) as executor:
        snv_max_pro = list(tqdm(executor.map(verify_snv_site, query_list), total=len(query_list)))

print(len(snv_max_pro))

snv_max_pro = np.array(snv_max_pro)
fake_snv_sites = to_be_verified_sites[snv_max_pro<cond_pro]
del snv_max_pro
final_snv = [i for i in range(num_snv) if i not in fake_snv_sites]

print(len(final_snv))

###  update snv_sites file
f = open(file_snv,'w')
for i in final_snv:
    f.write(str(snv_sites[i]+1)+'\t')


### check the satisfaction of SNV sites's number
if len(final_snv)<smallest_snv:
    print("No enough SNV sites.\nexit")
    f.write("\nA small number of SNV sites indicates one haplotype\nexit")
    f.close()
    exit()


f.close()
if only_snv != '0':
    exit()

###  update some variables
Coun = Coun[final_snv]
seq_mat = seq_mat[:,final_snv]
snv_base = snv_base[final_snv]
major_base = major_base[final_snv]
Qstart = []
Qend = []
seq_flag = []
snv_range = np.array([i for i in range(len(final_snv))])
for i in range(len(seq_mat)):
    temp1 = seq_mat[i]!='-'
    temp2 = np.logical_or(seq_mat[i]==snv_base,seq_mat[i]==major_base)
    if np.count_nonzero(temp1)==0:
        Qstart.append(snv_range[-1])
        Qend.append(snv_range[-1])
    else:
        Qstart.append(snv_range[temp1][0])
        Qend.append(snv_range[temp1][-1])
    seq_flag.append(temp2)

Qstart = np.array(Qstart)
Qend = np.array(Qend)
seq_flag = np.array(seq_flag)
index_r = np.argsort(Qstart)

print("reads graph construction")
file = f'{file_prefix}_reads_graph.txt'

def graph_construction(read_id,snv_overlap=ovlap_read, threshold=weight_read,Reads_all=Reads_all, Map_acc=Map_acc,Qstart=Qstart,Qend=Qend,seq_mat=seq_mat,seq_flag=seq_flag,index_r=index_r,m=len(index_r),step=sub_graph):
    index_i = index_r[read_id]
    read_name_i = Reads_all[index_i]
    qstart_i = Qstart[index_i]
    qend_i = Qend[index_i]
    qmid_i = int((qstart_i + qend_i)/2)
    seq_i = seq_mat[index_i]
    seq_flag_i = seq_flag[index_i]
    edge = []
    for j in range(read_id+step,m,step):
        index_j = index_r[j]
        qend_j = Qend[index_j]
        qstart_j = Qstart[index_j]
        read_name_j = Reads_all[index_j]
        if qstart_j > qmid_i:
            break
        seq_j = seq_mat[index_j]
        seq_flag_j = seq_flag[index_j]
        seq_flag_ij = np.logical_and(seq_flag_i, seq_flag_j)
        diff_flag = (seq_i[seq_flag_ij]!=seq_j[seq_flag_ij])
        d_score = np.count_nonzero(diff_flag)
        N = np.count_nonzero(seq_flag_ij)
        p_cor_i = Map_acc[read_name_i]
        p_cor_j = Map_acc[read_name_j]
        if N >= snv_overlap:
            weight_temp = np.round((p_cor_i*p_cor_j)*(1-(0.1+d_score)/(N+0.1)),4)
            if weight_temp >= threshold:
                edge.append(f"{index_i}\t{index_j}\t{weight_temp}\n")
    return edge


def writeedge(index_and_file):
    edge = graph_construction(index_and_file[0])
    with open(index_and_file[1],'a') as f:
        for e in edge:
            f.write(e)

if sub_graph == 1:
    file = f'{file_prefix}_reads_graph.txt'
    query_list = [[i,file] for i in range(len(index_r)-1)]
    with concurrent.futures.ProcessPoolExecutor(thread) as executor:
        results = list(tqdm(executor.map(writeedge, query_list), total=len(query_list)))
else:
    for k in range(sub_graph):
        file_index = k
        file = f'{file_prefix}{k}_reads_graph.txt'
        query_list = [[i+k,file] for i in range(0,len(index_r)-1,sub_graph) if (i+k)<len(index_r)-1]
 #       if k==sub_graph-1:
 #           temp = query_list[-1][0]+1
 #           for i in range(temp,len(index_r)-1):
 #               query_list.append([i,file])
        with concurrent.futures.ProcessPoolExecutor(thread) as executor:
            results = list(tqdm(executor.map(writeedge, query_list), total=len(query_list)))


data = {'seq_mat':seq_mat,'reads':Reads_all,'Qstart':Qstart,'Qend':Qend,'Coun':Coun,'snv_base':snv_base,'major_base':major_base}
f_pic = open(file_prefix+'_matrix.pickle','wb')
pickle.dump(data,f_pic,protocol=pickle.HIGHEST_PROTOCOL)
f_pic.close()

exit()
