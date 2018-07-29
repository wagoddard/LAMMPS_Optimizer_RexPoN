
import prmtop, bgf, amberFF, forceField, mol2, shortpar

ff = forceField.ForceField()
shortpar.read(ff)



'''
#read prm inp
AL = prmtop.readTopoCord('azobenzene.prmtop','azobenzene.inpcrd')
# make FF
FF = forceField.ForceField()
#read amber FF
amberFF.readAmberDat('/home/vcvicek/vasek/amber/gaff.dat',FF)
amberFF.readAmberFrcmod('azobenzene.frcmod',FF)
#write prm
prmtop.writeTopology('myazo.prmtop',AL,FF)
#write crd
prmtop.writeCoordinates('myazo.inpcrd',AL)
'''
'''
# make FF
FF = forceField.ForceField()
#read amber FF
amberFF.readAmberFF(FF,'/home/vcvicek/vasek/amber/gaff.dat',debug=True)
amberFF.readAmberFrcmod(FF,'atoms.frcmod',debug=True)

list = ['atoms1']
for a in list:
    #read from bgf
    AL = bgf.readFile(a+'.bgf')
    #write prm
    prmtop.writeTopology(a+'.prmtop',AL,FF,debug=True)
    #write crd
    prmtop.writeCoordinates(a+'.inpcrd',AL)
'''

