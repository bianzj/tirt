from old.physical.endmember_hom_hom import *
from old.physical.endmember_crown_hom import *
from old.variable import *


'''
树冠+树干+树下植被+土壤
其他场景：树冠+土壤
'''

variables = sample()
end_hom = Endmember_hom_hom()
end_hom.set_angle(variables)
end_hom.set_optical(variables)
end_hom.set_structure(variables)
end_hom.set_thermal(variables)
p_hom = end_hom.run()

plt.plot(p_hom,'-')
plt.show()

