import math as m

def ConservationScore(indata,alphabet,weight='probability',logbase=m.e): #calculates the information in a sequence of data made up from elements in alphabet. in our case, this is the information in a column of the alignment
    if weight=='probability':
        weight=1./float(len(alphabet))
    else:
        weight=weight
    counts=[]
    for i in alphabet:
        counts.append(sum([x==i for x in indata]))
    pweight=float(sum(counts))
    infos=[]
    for z in counts:
        if z==0:
            infos.append(0.)
        else:
            infos.append(float(z)/pweight*m.log(float(z)/pweight/weight,logbase))
    return sum(infos)