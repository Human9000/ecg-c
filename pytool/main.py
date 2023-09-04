import re
from cstruct import c_struct


dim0_vals = []

print("// ====== static array  refactor =========== ")
for l in c_struct.split('\n'):
    p = 0
    for i, v in enumerate(l):
        if v not in [' ', '\t']:
            break
        p = i + 1
    l = l[p:]
    r = l[:]
    l = re.sub('//.*', '', l)
    l = re.sub('static ', '', l)
    l = re.sub(';', '', l)
    l = re.sub('=.*', '', l)

    if len(l) == 0:
        # print()
        continue
    tva = l
    dim = tva.count("[")
    if dim == 0:
        tv = tva
    elif dim == 1:
        index = tva.index('[')
        tv = tva[:index]
    elif dim == 2:
        index = tva.index('[')
        tv, a = tva[:index], tva[index + 1:]
        index = a.index('[')
        a = a[index:]

    for i in range(len(tv)-1,-1,-1):
        if tv[i] == ' ':
            t, v = tv[:i], tv[i+1:]
            break

    if dim == 0:
        dim0_vals.append([t,v,r])
        # print(f"{t:10s} {v};")

    elif dim == 1:
        print(f"{t} *{v} = ctx->{v}; // {r}")
    elif dim == 2:
        print(f"{t} (*{v}){a} = ctx->{v}; //{r}")

print("\n// ====== static not array refactor =========== ")
for t,v,r in dim0_vals:
    print(f"{t:10s} {v}")

print("\n// ====== static refactor END =========== ")
正则匹配_添加访问符号 = [
    '([^*])(ctx_[a-zA-Z_0-9]+)',
    '$1(*$2)'
]
