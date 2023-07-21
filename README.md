# Systematics

When first setting up, or any time anything in Variation.* or CompositeDatat.* is changed, run these commands with a ROOT 6 installation (likely will not work on RACF though for reasons apart from ROOT):

```bash
./makeSharedLib.sh Variation
./makeSharedLib.sh CompositeData
```

Go into `calculateSystematics.cxx` and make sure the first argument for each `Variation` object is the path and prefix of your results variations. The normal results are also considered a 'variation' as well. To run the script and calculate systematic uncertainties, use 

```bash
root -l -b -q calculateSystematics.cxx
```

The `Variation` objects are then compined into one `CompositeData` object for each type of cut. Each `CompositeData` object can accept up to 5 `Variation` objects to combine. `CompositeData` does the background work to calculate the standard deviations and the significance of each cuts's systematic uncertainty contribution. All you need to do is use the `CompositeData` objects to pull out the significance value and the systematic uncertainty contribution.

Each `CompositeData` holds vectors corresponding to each important TProfile, for example, proton flow measurements vs centrality is called `p_vn_pp` (`p_` for 'profile', `v_` for 'vector'). To get the corresponding data that was compiled from all variations, use `compositeDataObject->v_vn_pp`. These vectors hold another custom data type called `DataPoint` defined within `CompositeData.h` which holds all necessary systematic info for each point.

## Barlow Check

The variable for every `DataPoint` that tells you the significance of the systematic contribution is called `deltaByDeltaError`. To determine if a particular cut and its variations (one `compositeDataObject`) produced a systematic uncertainty contribution that was significant at one bin of the proton flow vs centrality measurements, the command would be:

```c++
compositeDataObject->v_vn_pr.at(binNumber).deltaByDeltaError > 1.0
```

## Systematic Uncertainty Values

The systematic uncertainty contribution can be extracted from each `DataPoint` using `stdDev`, or the square of that with `variance`:

```c++
compositeDataObject->v_vn_pr.at(ithBin).stdDev
compositeDataObject->v_vn_pr.at(ithBin).variance
```

## Adding More Histograms

If there is a particular plot that you need to create systematic uncertainties for that is not present in `Variation` and `CompositeData`:

1) In `Variation.h`, create an instance of a `TH1D` corresponding to that plot.
2) In `Variation.cxx`, create a `delete` command in the destructor for that `TH1D` made in the previous step.
3) In `Variation.cxx` within the `initialize()` function, pull the histogram/profile from the file and convert it into a `TH1D` if it isn't already.
4) In `CompositeData.h`, create a `std::vector<DataPoint>` object corresponding to the histogram prepared in the last step.
5) In `CompositeData.h`, create an instance of a `TH1D` to hold the significance values related to the Barlow check.
6) In `CompositeData.cxx`, create a `delete` command in the destructor for that `TH1D` made in the previous step.
7) In `CompositeData.cxx` within the `initialize()` function, finish the construction of the `TH1D` made in the previous step. It should have identical bins and bounds as the flow plot it is related to.
8) In `CompositeData.cxx` within the `saveDetails()` function, add a call to `addRawValuesToFile()` and `addBarlowValuesToFile()` for the new plot you're trying to make. "detailsFile" should be the first argument to both and the vector made in step 4 should be the third argument to both. The second argument of `addRawValuesToFile()` is the histogram name from the normal version of the plot, and the second argument of `addBarlowValuesToFile()` is the histogram created in step 5 and 7.
9) In `CompositeData.cxx` within the `combine()` functions, add a call to `mergePoints()` for the new vector you made in step 4, where you supply that vector as the last argument. The first argument is the associated histogram from the normal data, and all of the middle arguments are the same histogram from the different variations. There are multiple versions of `combine()` and `mergePoints()` to facilitate cuts that have different numbers of variations available.
10) Run the commands mentioned at the top of this README file and everything should be compiled and ready to go. You can then retrieve the systematics for your new plot within `calculateSystematics.cxx` in the same way as the other plots.
