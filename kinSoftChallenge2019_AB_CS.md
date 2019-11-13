# kinSoftChallenge2019 Analysis
Anders Barth, Claus Seidel
## Analysis workflow
### Determination of FRET efficiencies
The minimal number of states and their FRET efficiencies are determined by Gaussian fitting of frame-wise FRET efficiency histograms.
### Analysis of kinetics
We employ fluorescence correlation spectroscopy to analyze the time traces of donor and fluorescence signal. Three different approaches are used:

1. Traditional "color-FCS", i.e. correlation of the raw intensity traces
2. "filtered-FCS" analysis - By preceding the analysis by a step-finding algorithm, states are identified and the resulting state trajectories are correlated

### Step-finding algorithm
We apply the algorithm presented by (Aggarwal, 2012)[][#Salapaka] to identify steps in the FRET efficiency traces. The algorithm does not assume any particular kinetic scheme but estimates the optimal number of steps based on the noise of the signal. Overfitting avoided by introducing a penalty for each transition. For the analysis, we set an estimated noise of $\sigma_E$ = 0.05. An exemplary result of the step-finding is shown in Figure 1.

/Test_Datasets/step_finding_example.png "Example of the step-finding algorithm. The idealized signal trajectory (black) is estimated from the noisy data (gray). Thresholds for the digitizations of the FRET efficiency trajectory into states are shown as dashed lines."

## FCS model function
### 1. Color-FCS
For two-state dynamics, analytical functions are known for the color-FCS curves.
### 2. step-finding FCS
In filtered-FCS, the correlation functions are obtained directly from the matrix exponential of the transition rate matrix:

## Analysis results
### Test datasets
#### FRET efficiency histograms
/Test_Datasets/FRET_combined.png (FRET efficiency histograms of the test datasets. Top row: Framewise FRET efficiency histograms and Gaussian fits. Weighted residuals are shown above. Bottom row: After step-finding, the FRET efficiency distribution narrows. Dashed lines indicated the applied thresholds to digitize the FRET efficiency trajectories for the calculation of state-filtered FCS curves.)

[#Salapaka]: T. Aggarwal, D. Materassi, R. Davison, T. Hays, M. Salapaka, Detection of Steps in Single Molecule Data. Cel. Mol. Bioeng. 5, 14–31 (2012).


