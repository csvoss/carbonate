"""
reaction_functions
Database-friendly reactions.
"""


import reactions as rxns
from toMolecule import moleculify as molec

magicify_it = lambda x: molec("[AuH8]")

hydrogenate_it = rxns.hydrogenate
radicalhydrobrominate_it = lambda x: rxns.radicalhydrohalogenate(x, "Br")

hydrobrominate_it_once = lambda x: rxns.hydrohalogenate1eq(x, "Br")
hydroiodinate_it_once = lambda x: rxns.hydrohalogenate1eq(x, "I")
hydrochlorinate_it_once = lambda x: rxns.hydrohalogenate1eq(x, "Cl")
hydrobrominate_it = lambda x: rxns.hydrohalogenate(x, "Br")
hydroiodinate_it = lambda x: rxns.hydrohalogenate(x, "I")
hydrochlorinate_it = lambda x: rxns.hydrohalogenate(x, "Cl")
brominate_it_once = lambda x: rxns.halogenate1eq(x, "Br")
iodinate_it_once = lambda x: rxns.halogenate1eq(x, "I")
chlorinate_it_once = lambda x: rxns.halogenate1eq(x, "Cl")
brominate_it = lambda x: rxns.halogenate(x, "Br")
iodinate_it = lambda x: rxns.halogenate(x, "I")
chlorinate_it = lambda x: rxns.halogenate(x, "Cl")

epoxidate_it = rxns.epoxidate

acidhydrate_it = lambda x: rxns.acidhydrate(x, molec("O"), False)
acidhydrate_it_hgso4 = lambda x: rxns.acidhydrate(x, molec("O"), True)
acidhydrate_it_ethanol = lambda x: rxns.acidhydrate(x, molec("CCO"), False)
acidhydrate_it_hgso4_ethanol = lambda x: rxns.acidhydrate(x, molec("CCO"), True)
acidhydrate_it_auto = lambda x: rxns.acidhydrate(x, x, False)
acidhydrate_it_hgso4_auto = lambda x: rxns.acidhydrate(x, x, True)

bromohydrate_it_water = lambda x: rxns.halohydrate(x, molec("O"), "Br")
bromohydrate_it_ethanol = lambda x: rxns.halohydrate(x, molec("CCO"), "Br")
bromohydrate_it_auto = lambda x: rxns.halohydrate(x, x, "Br")
iodohydrate_it_water = lambda x: rxns.halohydrate(x, molec("O"), "I")
iodohydrate_it_ethanol = lambda x: rxns.halohydrate(x, molec("CCO"), "I")
iodohydrate_it_auto = lambda x: rxns.halohydrate(x, x, "I")
chlorohydrate_it_water = lambda x: rxns.halohydrate(x, molec("O"), "Cl")
chlorohydrate_it_ethanol = lambda x: rxns.halohydrate(x, molec("CCO"), "Cl")
chlorohydrate_it_auto = lambda x: rxns.halohydrate(x, x, "Cl")

hydroborate_oxidate_it = rxns.hydroborate
hydroborate_oxidate_it_1 = rxns.hydroborate1
hydroborate_oxidate_it_2 = rxns.hydroborate2

dihydroxylate_it = rxns.dihydroxylate
ozonolyse_it = rxns.ozonolyse
sodium_ammonia_it = rxns.sodiumAmmonia
lindlar_it = rxns.lindlar
alkyne_deprotonate_it = rxns.alkyneDeprotonate
tert_butoxide_it = rxns.tertButoxide
acetylide_add_it = lambda x: rxns.acetylideAdd(x, x)
