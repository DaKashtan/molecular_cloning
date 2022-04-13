from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

vst_numb = 3
vstavki = []
for i in range(3):
    vst = open(input('file:'), 'r')
    vst_seq = parse(vst, 'fasta')
    for j in vst_seq:
        print(j)
        vst1 = j.seq
        v = str(vst1) + '!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    vstavki.append(v)
    vstavki.append('>>>>>')
    vst.close()
print(vstavki)