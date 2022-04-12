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

#сюда функцию праймер для варианта 1

def primer(n):
    primers=[]
    primersL = good_sitesL[n] + vstavki_sort[n][0:21]
    primers.append(primersL)
    primersRb = ''.join(reverse_complementary(vstavki_sort[n][-1:-16:-1]))
    primersR = good_sitesL[n] + good_sitesL[n+1] + primersRb #к обратному праймеру добавляем первый + второй сайты
    primers.append(primersR)
    return primers

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
print(good_sitesL)

#проверяем наличие участка множественного клонирования MCS
sites_distance = []
for i in good_sitesL:
    sites_distance.append(vector_sequence.index(i))
print(sites_distance)

MCS = True
for i in sorted(sites_distance):
    col = len(vstavki_sort)
    c = 0
    if abs(i - (i+1)) < 100:
        c +=1
if c < col:
    MCS = False
print(MCS)

if MCS:
    numb_restr = len(vstavki_sort)
    i=0
    result=vector_sequence
    while i!=numb_restr:
        site_restr=good_sitesL[i]
        a, b = result.split(site_restr,1)
        result = a+site_restr+vstavki_sort[i]+site_restr + b #теперь берем сайт рестриции и добавляем вставку
        i+=1
    #здесь код с праймерами
    print(result)

#второй вариант
else:
    for i in vstavki_sort:
        print(primer(vstavki_sort.index(i))) #списки праймеров для варианта 2 для каждой вставки


    ###### первая вставка = сайт рестрикции 1 + сама вставка + сайт рестрикции 2 + сайт рестрикции1
    unic_sites = []
    for i in sites:
        if i not in vector_sequence:
            unic_sites.append(i)  #список сайтов, не встречающихся в векторе
            if i in vstavki:
                unic_sites.remove(i)
    print(unic_sites)

    numb_restr = len(vstavki_sort)
    a, b = vector_sequence.split(good_sitesL[0], 1)
    result = a + good_sitesL[0] + vstavki_sort[0] + unic_sites[0] + good_sitesL[0] + b #вставили первую вставку, добавив за ней уникальный сайт

    i=0
    while i!=numb_restr-1:
        a, b = result.split(unic_sites[i],1)
        result = a+unic_sites[i]+vstavki_sort[i+1]+ unic_sites[i+1] + unic_sites[i] + b #теперь берем сайт рестриции и добавляем вставку с сайтами
        i+=1
print(result)




# вывод манипуляций
#new_file = input('Введите название нового файла ')
#my_file = open(new_file, "w")
#a='ddd'
#print(my_file.write("Закажите праймер "+ a +'\n') ) #вставить переменную праймера
#print(my_file.write("Возьмите рестриктазу " + a +'\n')) #вставить переменную из первого выхода цикла, который ищет сайты
#print(my_file.write("Порежьте вектор рестриктазой " + a +'\n'))#вставить переменную из первого выхода цикла, который ищет сайты
#print(my_file.write("Возьмите рестриктазу " + a +'\n') )
#print(my_file.write('Возьмите порезанный вектор и вставку ' +a +'\n')) #переменная из части с праймерами.
#print(my_file.write('Лигируйте ' + a +'\n'))
#print(my_file.write('Проведите рестрикционный анализ.' + a +'\n'))
#my_file.close()


vector.close()
vst.close()
sites1.close()



