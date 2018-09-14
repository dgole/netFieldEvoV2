#!/bin/bash

module purge
module load gcc
module load python/3.6.1

# means still needs to be run
## means done

# 
#python3 ../python/netFieldEvoSSExors.py 30110 
#python3 ../python/netFieldEvoSSExors.py 30111
#python3 ../python/netFieldEvoSSExors.py 30112 
#python3 ../python/netFieldEvoSSExors.py 30113
#python3 ../python/netFieldEvoSSExors.py 30114

#python3 ../python/netFieldEvoSSExors.py 30120 
#python3 ../python/netFieldEvoSSExors.py 30121
#python3 ../python/netFieldEvoSSExors.py 30122 
#python3 ../python/netFieldEvoSSExors.py 30123
#python3 ../python/netFieldEvoSSExors.py 30124

#python3 ../python/netFieldEvoSSExors.py 30130 
#python3 ../python/netFieldEvoSSExors.py 30131
#python3 ../python/netFieldEvoSSExors.py 30132 
#python3 ../python/netFieldEvoSSExors.py 30133
#python3 ../python/netFieldEvoSSExors.py 30134

#python3 ../python/netFieldEvoSSExors.py 30200
#python3 ../python/netFieldEvoSSExors.py 30201

#python3 ../python/netFieldEvoSSExors.py 30210
#python3 ../python/netFieldEvoSSExors.py 30211

#python3 ../python/netFieldEvoSSExors.py 30220
#python3 ../python/netFieldEvoSSExors.py 30221

#python3 ../python/netFieldEvoSSExors.py 30231


#python3 ../python/netFieldEvoSSExors.py 30230


# DZ 160

#python3 ../python/netFieldEvoSSExors.py 30400 # -1.e3, +1.e-3
#python3 ../python/netFieldEvoSSExors.py 30401 # +1.e3, -1.e-3
#python3 ../python/netFieldEvoSSExors.py 30402 # -1.e3, +1.e-2
#python3 ../python/netFieldEvoSSExors.py 30403 # +1.e3, -1.e-2
#python3 ../python/netFieldEvoSSExors.py 30404 # -1.e3, +1.e-1
#python3 ../python/netFieldEvoSSExors.py 30405 # +1.e3, -1.e-1

# DZ 160, tFlip=10
#python3 ../python/netFieldEvoSSExors.py 31000  # init 1.e-3  # looks pretty bad
#python3 ../python/netFieldEvoSSExors.py 31001  # init 1.e-2  # looks ok, maybe run this suuuper long like 100 or 200 cycles

python3 ../python/netFieldEvoSSExors.py 31002  # init 1.e-2  # longer

# DZ 0
#python3 ../python/netFieldEvoSSExors.py 31010  # tFlip=3   # looks meh, maybe run longer?
#python3 ../python/netFieldEvoSSExors.py 31011  # tFlip=10  # looks really good

# ramp up
#python3 ../python/netFieldEvoSSExors.py 31020 # +1.e3, +1.e-2

#python3 ../python/netFieldEvoSSExors.py 31021 # ramp up, B -1.e3 -1.e-6, tflip 1
python3 ../python/netFieldEvoSSExors.py 31022 # -1.e3, -1.e-6 tflip 10

# DZ 200
python3 ../python/netFieldEvoSSExors.py 31020  # init 1.e-2  



# DZ 200 31030
# DZ 220 31031
# DZ 240 31032
# DZ 160 31033
# DZ 120 31034
# DZ 120 31035 3 years

















