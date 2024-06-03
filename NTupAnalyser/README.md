# NTupAnalyser

This code will analyse the NTuples which were obtained by the HGam package. It will set the right c-tagging and c-tagging efficiencies also will make the histograms and save them in the file.

Also, it will do the optimization studies and will save all the nessesary hiatograms, which will be used by the plotting functions

## Quik start:

**To clone the git repository**

```
git clone ssh://git@gitlab.cern.ch:7999/aivina/NTupAnalyser.git

```

**Get the latest changes**:

```
git pull --all

```

**To check what you have done after you made changes**:
```
git status
```

**To add changes**:

```
git add filename.txt

```

**Commit the changes**:
```
git commit -m "a short message here"
```

**Send the changes to the brunch**:
```
git push -u origin master
```
```
Now your changes in the git repository and you need to merge them with the master branch:

git checkout master
git merge my-local-branch
git push origin master
```

## How to run the code:

asetup AnalysisBase,21.2.40,here

cmake ../source   and then cmake --build .

source $TestArea/../build/$AnalysisBase_PLATFORM/setup.sh

To run the code: run_ntuple --inputDir  --Nevents 

**The options for running**

```
-- output_file - the name of your output file

-- inputDir - the path to the directory where your inputs are stored

-- submitDir - the name of your output folder 

-- histname - the name of the CutFlow histogram in the MxAOD

-- RootFile - the namde of the Root file which the output of the HGam

-- mcname - the name of the Mc sample which we are dealing with (useful for the creating output txt file)

-- Nevents - Maximum number of events

-- reco_jet_collection - Jet collection which you prefer (EMTopo or PFlow)

-- ispTcut - dO you want to apply pT cut on photons? yes -> true

-- phPtcut - which cut to set on the photons

-- doCutFlow - do you want to make the cutflow histogram? say yes -> true

-- doOptimization - the flag which will run over jets to optimize the ctagging point

-- bfrac - will be used latr, when the optimization is done

-- jetpT - what the pT of jet you want to set?

-- jeteta - what eta of jet you want to set?

-- bfrac - this is the mistag efficiency for c tagging

-- isSignal - so far I am keeping it (remove when fix the Cross section in the McSamples.config)

-- DL1cut - the cut to make the c tagging

-- btagger -  what b tagger to use? MV2c10 or DL1 with different eff cut, to choose please do MV2c10_FixedCutBeff_* or DL1_FixedCutBeff_*

-- InDS - the name of the inputs on the Grid

-- OutDS - the name of the outputs from the Grid

-- isGridJob - is it a Grid job or not ? Default 0

-- nFilesPerJob - Number of files per job to be submitted on the Grid

```