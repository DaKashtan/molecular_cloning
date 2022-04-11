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
        if vector_sequence.count(i) == 1:
            good_sites[i] = vector_sequence.count(i)  #ищем все сайты с уникальным вхождением в вектор
print(good_sites)
good_sitesL = list(good_sites.keys())
good_site = good_sitesL[0] #берем первый из них

def primer(n):
    primers=[]
    primersL = good_sitesL[n] + vstavki_sort[n][0:21]
    primers.append(primersL)
    primersRb = ''.join(reverse_complementary(vstavki_sort[n][-1:-16:-1]))
    primersR = good_sitesL[n] + good_sitesL[n+1] + primersRb #к обратному праймеру добавляем первый + второй сайты
    primers.append(primersR)
    return primers
for i in vstavki_sort:
    print(primer(vstavki_sort.index(i)))


###### первая вставка = сайт рестрикции 1 + сама вставка + сайт рестрикции 2 + сайт рестрикции1
vst1 = good_site + vstavki_sort[0] + good_sitesL[1] + good_site
print(vst1)

def vstavka_form(n):
    vstavka = good_sitesL[n] + vstavki_sort[n] + good_sitesL[n+1]
    return vstavka

#### режем (вены)
#part_vector_1, part_vector_2 = vector_sequence.split(good_site) #режем вектор
#res1 = part_vector_1 + vst1 + part_vector_2 #вставляем туда вставку
#print(res1)

### теперь вторая вставка - функция с праймерами и короче всё то же самое, только грамотно зациклить надо НО! с каждой
#новой итерацией нужно проверять, нет ли сайта в предыдущей вставке




#for v in vstavki_sort:
    #vst_for_insertion = {}
    #vst_for_insertion[vstavki_sort.index(v)] = primer(v)
    #part_vector_1, part_vector_2 = vector_sequence.split(good_sitesL[n])
    #res1 = part_vector_1 + vst_for_insertion[n] + part_vector_2
#print(vstavka_form(vstavki_sort[0]))

