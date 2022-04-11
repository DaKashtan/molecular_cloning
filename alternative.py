from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

def reverse_complementary(a):
    comp=[]
    for i in a:
        if i =='A':
            comp.append('T')
        elif i =='T':
            comp.append('A')
        elif i =='G':
            comp.append('C')
        elif i =='C':
            comp.append('G')
    comp = comp[-1::-1]
    return comp

#vector = open(input("Задайте адрес к файлу с векторной последовательностью: "),'r') #последовательность вектора
vector = open("C:/Users/yeba/Desktop/vector.fasta", 'r')
records = parse(vector, "fasta")
for record in records:
    vector_sequence = record.seq
vector_sequence=str(vector_sequence)
print(vector_sequence)


vstavki = []
#vst = open(input("Задайте адрес к файлу с вставками: "), 'r') #открываем файл со вставками
vst = open("C:/Users/yeba/Desktop/vstavki.txt", 'r')
vstavki = [i for i in vst.read().splitlines() if i]
order = list(map(int, input('Введите порядок организации вставок в векторе: ').split()))
count_order=1
numb_order=len(order)
vstavki_sort=[]
while count_order!=len(order)+1:
    ind_order=order.index(count_order)
    vstavki_sort.append(vstavki[ind_order])
    count_order+=1
print(vstavki_sort)


#открываем файл с сайтами рестрикции
#sites1 = open(input("Задайте адрес к файлу с сайтами рестрикции: "), 'r') ## правка
sites1 = open("C:/Users/yeba/Desktop/resrtr.txt", 'r')
sites = [i for i in sites1.read().splitlines() if i] #получаем список с сайтами
print(sites)

good_sites= {}
for i in sites:
    if i in vector_sequence:
        if vector_sequence.count(i) <2:
            good_sites[i] = vector_sequence.count(i)  #ищем все сайты с уникальным вхождением в вектор
print(good_sites)
good_sitesL = list(good_sites.keys())
good_site = good_sitesL[0] #берем первый из них

###получаем праймеры к первой вставке ВОТ ЭТО СДЕЛАТЬ В ЦИКЛ ДЛЯ ВСЕХ ВСТАВОК (а лучше функцию написать)
primers1L = good_site + vstavki_sort[0][0:21]
primers1Rb = ''.join(reverse_complementary(vstavki_sort[0][-1:-16:-1]))
primers1R = good_site + good_sitesL[1] + primers1Rb #к обратному праймеру добавляем первый + второй сайты
print(primers1L)
print(primers1R)

###### первая вставка = сайт рестрикции 1 + сама вставка + сайт рестрикции 2 + сайт рестрикции1
vst1 = good_site + vstavki_sort[0] + good_sitesL[1] + good_site
print(vst1)


