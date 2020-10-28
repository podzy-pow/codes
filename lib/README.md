# build
build requires `pybind11` version >= 2.5.  
```
 python setup.py build_ext --inplace
 ```
 
 # Galois field(GF) 4
The library operates with codes in GF4, denoting w=u and w̄=v of traditional notation. The tables of operations in the field are
<table>
<tr><td>

|x|0|1|w|w̄|
|-|-|-|-|-|
|0|0|0|0|0|
|1|0|1|w|w̄|
|w|0|w|w̄|1|
|w̄|0|w̄|1|w|

</td><td>

|+|0|1|w|w̄|
|-|-|-|-|-|
|0|0|1|w|w̄|
|1|1|0|w̄|w|
|w|w|w̄|0|1|
|w̄|w̄|w|1|0|

</td></tr> </table>

 # Usage example

 ```python
>> from codeslib import *
>> gf4('v') * gf4('u')
1
>> c = Code([Series("0u"), Series("0v"), Series("11")], 3, 1)
>> c.isSelfOrthogonal()
False
>> c.minDistance()
4
>> c = Code([Series("11"), Series("1u"), Series("1v")], 3, 1)
>> c.isSelfOrthogonal()
True
>> c.minDistance()
6
>> c.findOrthogonal()
00|1v|uv||1|v|u
 ```

1/3-code search with degree of 6

 ```python
from codeslib import *
s = SearchSelfOrthogonal(3, 6)
codes = s.find()
codes.sort(key=lambda c: c.weight())
print('Codes found, searching for optimal ones')
print('Total number of codes:', len(codes))
vlist = []
for i, code in enumerate(codes):
    orth_code = code.findOrthogonal()
    u, v = orth_code.minDistance(), code.minDistance()
    if u != 0:
        vlist.append((u, v))

miu, mau = min(x for x, y in vlist), max(x for x, y in vlist)
miv, mav = min(y for x, y in vlist), max(y for x, y in vlist)

table = [[0] * (mav-miv+1) for _ in range(miu, mau+1)]

for x, y in vlist:
    table[x-miu][y-miv] += 1

s = ' ' * 5 + '||' + ' | '.join('{:5}'.format(x) for x in range(miv, mav+1))
print(s)
print('-' * len(s))
for i, row in enumerate(table):
    print('{:5}||'.format(miu + i) + ' | '.join('{:5}'.format(x) for x in row))
``` 
