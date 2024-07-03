import sys

class Transcript():
    def __init__(self, trans_id, chrom, start, end, anntype, strand):
        self.trans_id = trans_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.anntype = anntype
        self.strand = strand


query_dict = {}
with open("targ_f.tsv", 'r') as f:
    """
    ENST_ID START END
    """
    for line in f:
        orig_id, enst_id, gstart, gend = line.strip().split()
        query_dict.setdefault(enst_id, []).append((int(gstart), int(gend), orig_id))

#print(query_dict)
trans_dict = {}
with open("/slurm/home/yrd/linlab/suxinwan/data/reference/ensembl/Homo_sapiens.GRCh38.111.chr.gtf", 'r') as f:
    for line in f:
        # 跳过以 "#" 开头的行
        if line.startswith('#'):
            continue
        tmp = line.strip().split("\t")
        chrom = tmp[0]
        anntype = tmp[2]
        if not anntype in ["exon"]:
            continue
        start = int(tmp[3]) 
        end = int(tmp[4])
        strand = tmp[6]
        info = tmp[8]
        for item in info.split(";"):
            if not item:
                break
            tmp2 = item.split()
            if tmp2[0] == "transcript_id":
                trans_id = tmp2[1].strip('"')
                if trans_id in query_dict:
                    trans_dict.setdefault(trans_id, []).append(Transcript(trans_id, chrom, start, end, anntype, strand))
                break

for trans_id, seg_list in trans_dict.items():
    seg_list.sort(key = lambda x:x.start)

"""
for ll in trans_dict.values():
    for t in ll:
        print(t.trans_id, t.chrom, t.start, t.end, t.strand, t.anntype)
    break
"""
print("ID", "CHROM", "GENOME_START", "GENOME_END", "TRANS_START", "TRANS_END", "STRAND", sep='\t')
with open("output.txt", "w") as f:
    f.write("ORF_ID\tTrans_ID\tCHROM\tGENOME_START\tGENOME_END\tTRANS_START\tTRANS_END\tSTRAND\n")
    for trans_id, pos_tuples in query_dict.items():
        for pos_tuple in pos_tuples:
            qstart, qend, orig_id = pos_tuple
            trans_seg_list = trans_dict[trans_id]
            trans_start = -1
            trans_end = -1
            pre_seg_len = 0
            for t in trans_seg_list:
                if t.start <= qstart <= t.end:
                    trans_start = qstart - t.start + 1 + pre_seg_len
                if t.start <= qend <= t.end:
                    trans_end = qend - t.start + 1 + pre_seg_len
                pre_seg_len += t.end - t.start + 1

            if trans_dict[trans_id][0].strand == "-":
                trans_start = pre_seg_len - trans_start + 1
                trans_end = pre_seg_len - trans_end + 1

            f.write(f"{orig_id}\t{trans_id}\t{t.chrom}\t{qstart}\t{qend}\t{trans_start}\t{trans_end}\t{t.strand}\n")