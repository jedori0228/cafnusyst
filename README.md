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

## Weights-only (friend-tree) output

By default `UpdateReweight` writes a full flattened CAF with the recomputed
weights attached, duplicating all of the unchanged content. With
`--weights-only` it instead writes a slim file whose `StandardRecord`s have
only `mc.nu[i].syst_dials` populated and everything else default-constructed.
ROOT compresses the repeated empty branches, so the file stays small while
keeping the standard CAF schema, and it is entry-aligned with the input CAF so
it can be used as a friend. The GENIE tree is not written in this mode; the
`globalTree` still is, so the weight indices remain interpretable.

```
UpdateReweight -c weighter.fcl -i input_cafs.txt -o weights_friend.root --weights-only
```

An event cap (`-N`) is rejected in this mode, since a truncated output cannot
stay entry-aligned with the full parent CAF.

How to read it back depends on the `StandardRecord` branch name, set with
`--sr-branch` (default `rec`). Same-named branches in a friend are shadowed by
the host tree, so a naive `AddFriend` on a `rec` branch would silently read the
parent's (empty) weights:

- `--sr-branch rec` (default): open both files and step one `StandardRecordProxy`
  per tree in lockstep by entry; do not `AddFriend`.
- `--sr-branch <other>`: `AddFriend` resolves cleanly by bare name; read the
  weights under that prefix.

## Output provenance

Every output file records how it was written as two top-level `TNamed`
objects, so a reader can pick the right access pattern without having to know
the weight updater command that was used:

- `cafnusyst_srbranch` — title is the output `StandardRecord` branch name
  (e.g. `rec`), i.e. the value of `--sr-branch`.
- `cafnusyst_mode` — title is `weights-only` or `full`.

Read them back with the templated `TDirectory::Get<T>()` (returns `nullptr`
on a type mismatch):
```cpp
TString branch = f->Get<TNamed>("cafnusyst_srbranch")->GetTitle();
TString mode   = f->Get<TNamed>("cafnusyst_mode")->GetTitle();
```



