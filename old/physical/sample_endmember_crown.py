from old.physical.endmember_crown_hom import *
from old.variable import *

'''
树冠+树干+树下植被+土壤
其他场景：树冠+土壤
'''

variables = sample()
end_crown = Endmember_crown_hom()
end_crown.set_angle(variables)
end_crown.set_optical(variables)
end_crown.set_structure(variables)
end_crown.set_thermal(variables)


rad = end_crown.run()
plt.plot(rad,'-o')
plt.show()

