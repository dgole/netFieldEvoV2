#!/bin/bash

module purge
module load gcc
module load python/3.6.1

# each line runs 16 integrations
# each output folder is about 500 MB -->  Gigs per line here



./runEachMdotTcycle.sh 10 0.1
./runEachMdotTcycle.sh 11 0.3
./runEachMdotTcycle.sh 12 0.5
./runEachMdotTcycle.sh 13 0.7
./runEachMdotTcycle.sh 14 1.0












#########################################################################################################################
# alphaActive = 0.03
# alphaDead   = 0.001
#########################################################################################################################
#                        runIdBase bStar  mdot0   tCycle nCycles
##./runEachStellarMass.sh  100       -1.e3  1.e-7   6      40  # no cycles
##./runEachStellarMass.sh  101       -1.e3  1.e-8   6      40  # no cycles
##./runEachStellarMass.sh  102       -1.e3  1.e-9   6      40  # no cycles
##./runEachStellarMass.sh  103       -1.e3  1.e-10  6      40  # no cycles

##./runEachStellarMass.sh  110       -1.e3  1.e-7   20     40  # no cycles
##./runEachStellarMass.sh  111       -1.e3  1.e-8   20     40  # no cycles
##./runEachStellarMass.sh  112       -1.e3  1.e-9   20     40  # no cycles
##./runEachStellarMass.sh  113       -1.e3  1.e-10  20     40  # tiniest of cycles on Mstar=0.1

##./runEachStellarMass.sh  120       -1.e3  1.e-7   60     40  # no cycles
##./runEachStellarMass.sh  121       -1.e3  1.e-8   60     40  # tiny cycles on Mstar=0.1
##./runEachStellarMass.sh  122       -1.e3  1.e-9   60     40  # tiny cycles on Mstar=0.3, 0.1
##./runEachStellarMass.sh  123       -1.e3  1.e-10  60     40  # tiny cycles on Mstar=0.7, weak 0.3, full 0.1

##./runEachStellarMass.sh  130       -1.e3  1.e-7   200    40  # really tiny on Mstar=0.3 and lower
##./runEachStellarMass.sh  131       -1.e3  1.e-8   200    40  # really tiny on Mstar=0.7, weak at 0.3, 0.1
##./runEachStellarMass.sh  132       -1.e3  1.e-9   200    40  # tiny at 1.0, 0.7, weak at 0.5, factor of 4 or 5 at 0.3, solid at 0.1
##./runEachStellarMass.sh  133       -1.e3  1.e-10  200    40  # some at all mass




































