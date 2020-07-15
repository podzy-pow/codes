# build
build requires `pybind11` version >= 2.5.  
```
 python setup.py build_ext --inplace
 ```

 # usage example

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
 ```
 