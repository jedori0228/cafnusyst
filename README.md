# cafnusyst

This package provides tools to produce new flattened CAF files with updated caf::SRTrueInteraction::wgt including new cross-section reweights evaluated by [nusystematics](https://github.com/NuSystematics/nusystematics) from original input CAF files.

# dependencies

## duneanaobj

`duneanaobj` version `v03_15_00` or later

```
setup duneanaobj v03_15_00 -qe26:prof
```

## nusystematics

### Installation

```
# Dependencies
setup cmake v3_27_4
setup genie v3_04_02 -qe26:prof
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup boost v1_82_0 -qe26:prof

setup eigen v23_08_01_66e8f
setup fhiclcpp v4_18_04 -qe26:prof
```

```
# #{mywd} is your working area
cd ${mywd} # go to your working directory
mkdir nusystematics; cd nusystematics
git clone git@github.com:NuSystematics/nusystematics.git nusystematics-src
mkdir build; cd build
cmake ../nusystematics-src/
make install
```

### setup script

Whenever you open a new shell, run
```
source ${mywd}/nusystematics/build/Linux/bin/setup.systematicstools.sh
source ${mywd}/nusystematics/build/Linux/bin/setup.nusystematics.sh
```

# build

```
# #{mywd} is your working area
cd ${mywd} # go to your working directory
mkdir cafnusyst; cd cafnusyst;
git clone git@github.com:jedori0228/cafnusyst.git cafnusyst-src
mkdir build; cd build
cmake ../cafnusyst-src/
make install
```

# setup script

Whenever you open a new shell, run
```
source ${mywd}/cafnusyst-src/build/Linux/bin/setup.cafnusyst.sh
```

# Running UpdateReweight

You need a nusyst configuration fhicl file. An example is located in `example/zexpansion_weighter.ParameterHeader.fcl`.

Then write a txt file that contains the list of input CAF files:
```
$ cat input_cafs.txt
path/to/your/input/caf.root
```

Then run `UpdateReweight` using following command:
```
UpdateReweight -c example/zexpansion_weighter.ParameterHeader.fcl -i input_cafs.txt -o output_flat.caf.root
```



