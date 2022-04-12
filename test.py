from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

vector = open(input("Задайте адрес к файлу с векторной последовательностью: "),'r') #последовательность вектора
records = parse(vector, "fasta")
for record in records:
    vector_sequence = record.seq
vector_sequence=str(vector_sequence)
print(vector_sequence)


vstavki = []
vst = open(input("Задайте адрес к файлу с вставками: "), 'r') #открываем файл со вставками
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
sites1 = open(input("Задайте адрес к файлу с сайтами рестрикции: "), 'r')
sites = [i for i in sites1.read().splitlines() if i]
print(sites)

# Подбираем сайты рестрикции и сортируем в порядке включения в последовательность для организации вставок в необходимом порядке
for i in sites:
    count = 0
    if i in vstavki:
            sites.remove(i)  # удаляем из списка сайты, которые есть в вставках
good_sites= {}
for i in sites:
    if i in vector_sequence:
        if vector_sequence.count(i) == 1:
            good_sites[i] = vector_sequence.find(i)  #ищем все сайты с уникальным вхождением в вектор

good_sitesL = sorted(good_sites,reverse=True)
numb_restr=len(vstavki_sort)
i=0
result=vector_sequence
while i!=numb_restr:
    site_restr=good_sitesL[i]
    a, b = result.split(site_restr,1)
    result = a+site_restr+vstavki_sort[i]+site_restr + b #теперь берем сайт рестриции и добавляем вставку
    i+=1

print(result)

# вывод манипуляций
print('1. Проведение ПЦР... сайт рестрикции', site_restr)
print('2. Обработка рестриктазой...')
print('3. Лигирование...')
print('4. Конечный результат...', result)


vector.close()
vst.close()
sites1.close()



