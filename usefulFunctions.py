def printArray(a):
    string=[]
    for row in a:
        string.append(str(row).replace('[',' ').replace(']',' '))
    return string


def sortList(seq,slist):
    return [slist[i] for (v,i) in sorted((v,i) for (i,v) in enumerate(seq))]
