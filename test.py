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
sites1 = open(input("Задайте адрес к файлу с сайтами рестрикции: "), 'r') ## правка
sites = [i for i in sites1.read().splitlines() if i] #получаем список с сайтами
print(sites)

good_sites= []
for i in sites:
    if i in vector_sequence:
        good_sites.append(i) #добавляем в список гуд_сайтс сайты, которые есть в векторе
for i in good_sites:
    count=0
    while count!=len(vstavki):
        if i in vstavki[count]:
            good_sites.pop(count) #удаляем из списка сайты, которые есть в вставках
        count+=1
print(good_sites)

# Ищем количество включений  сайтов в вектор
sum_inclusions=[]
for i in good_sites:
    inclusions=re.findall(i,  vector_sequence)
    sum_inclusions.append(len(inclusions)) #Подсчет количества включений "хороших" сайтов
    inclusions=[]
numb_restr=len(vstavki)
i=0
while i!=numb_restr:
    ind_site=sum_inclusions.index(min(sum_inclusions))
    site_restr=good_sites[ind_site] #Находим сайт рестрикции
    next_site=good_sites[sum_inclusions.index(min(sum_inclusions))]
    a, b = vector_sequence.split(site_restr,1)
    result = a+site_restr+vstavki_sort[i]+site_restr + b #теперь берем сайт рестриции и добавляем вставку
    sum_inclusions.pop(ind_site)
    good_sites.pop(ind_site)
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



