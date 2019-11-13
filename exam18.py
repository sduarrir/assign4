#e2

def int_divisor(num):
    divlist = []
    for i in range(num -1, 1, -1):
        if (num % i == 0):
            divlist.append(i)
    divlist.append(1)
    return divlist


def perfectint(num):
    div = int_divisor(num)
    sum = 0
    for i in div:
        sum += i
    return sum == num
a = perfectint(9)
print(a)

def fastasplit(fastafile):
    inputfile = open(fastafile, "rt")
    outputfile = open('divfile.fastq', "rt")
    while True:
        tag = inputfile.readline()
        if not tag: break  # End of file, there is no more lines # Sequence
        seq = inputfile.readline().strip()
        mid = len(seq)/2  #  longitud of the sequence
        # '+'
        plusline = inputfile.readline()  # + line
        # quality line
        qline = inputfile.readline().strip()  # qualities line remove end of line
        #file 1
        outputfile.write('{0}{1}\n{2}{3}\n'.format(tag, seq[:mid], plusline[:mid], qline))  # print sequence
        outputfile.write('{0}{1}\n{2}{3}\n'.format(tag, seq[mid:], plusline[mid:], qline))  # print sequence
    inputfile.close()
    outputfile.close()


def mergeseq (seq):
    prev = seq[0]
    newsew = list(seq[0])
    for el in seq[1:]:
        if el != prev:
            newsew.append(el)
            prev = el
    return ''.join(newsew)

def recseq(seq,prev=None):
    if len(seq) == 0:
        return ""
    if prev:
        if prev == seq[0]:
            return recseq(seq[1:], prev)
    return seq[0]+recseq(seq[1:], seq[0])


def merge(left, right):
    if left[-1] == right[0]:
        return left[:-1] + right
    return left + right

def divseq (seq):
    if len(seq) <= 1:
        return seq
    mid = len(seq)//2
    left, right = divseq(seq[:mid]), divseq(seq[mid:])
    return merge(left, right)



seq = 'AAAAAIIIKKKKJNGHLLLLL'
m = divseq(seq)

print(m)