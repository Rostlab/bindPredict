from ProfParser import ProfParser

file = 'Q9Y530.reprof'
file2 = 'Q9Y530.prof'
pp2 = ProfParser(file2, method='profphd')
pp = ProfParser(file, method='reprof')

print(pp.seq)
print(len(pp.seq))
print(len(pp.p_rel_solv))

print(pp2.seq)
print(pp2.p_rel_solv)
print(pp2.o_rel_solv_10)

print(pp.o_rel_solv_10)
