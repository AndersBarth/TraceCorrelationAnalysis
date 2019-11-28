# kinSoftChallenge2019 Analysis
Anders Barth, Claus Seidel

## Analysis workflow
All analysis was done in MATLAB. The code is hosted in a [git repository](https://github.com/AndersBarth/TraceCorrelationAnalysis).
### Determination of FRET efficiencies
The minimal number of states and their FRET efficiencies are determined by Gaussian fitting of frame-wise FRET efficiency histograms. We used a Gaussian mixture model as implemented in the `fitgmdist`, based on an iterative Expectation-Maximization algorithm of the likelihood function.
### Analysis of kinetics
We employ fluorescence correlation spectroscopy to analyze the time traces of donor and fluorescence signal. Three different approaches are used:

1. Traditional "color-FCS", i.e. correlation of the raw intensity traces
2. "filtered-FCS" analysis - By preceding the analysis by a step-finding algorithm, states are identified and the resulting state trajectories are correlated

### Step-finding algorithm
We apply the algorithm presented by (Aggarwal, 2012) [1] to identify steps in the FRET efficiency traces. The algorithm does not assume any particular kinetic scheme but estimates the optimal number of steps based on the noise of the signal. Overfitting avoided by introducing a penalty for each transition. For the analysis, we set an estimated noise of based on the distribution width obtained from the Gaussian fitting analysis ($\sigma_E$ = 0.05-0.1). An exemplary result of the step-finding is shown in Figure 1.

![Figure 1: Example of the step-finding algorithm. The idealized signal trajectory (black) is estimated from the noisy data (gray). Thresholds for the digitizations of the FRET efficiency trajectory into states are shown as dashed lines.](Test_Datasets/step_finding_example.png)

After the step finding, the stepwise FRET efficiency histograms were examined to identify thresholds to digitize the FRET efficiency trajectory.

![Figure 2: FRET efficiency histograms of challenge dataset 2 before (left) and after (right) applying the step finding. Dashed lines indicate the used thresholds.](Challenge%20Datasets/sim_level2_final_publish_E_combined.png)

## Analysis results
### Test datasets
#### FRET efficiency histograms

![FRET efficiency histograms of the test datasets. Top row: Framewise FRET efficiency histograms and Gaussian fits. Weighted residuals are shown above. Bottom row: After step-finding, the FRET efficiency distribution narrows. Dashed lines indicated the applied thresholds to digitize the FRET efficiency trajectories for the calculation of state-filtered FCS curves.](Test_Datasets/FRET_combined.png)

[1]: T. Aggarwal, D. Materassi, R. Davison, T. Hays, M. Salapaka, Detection of Steps in Single Molecule Data. Cel. Mol. Bioeng. 5, 14–31 (2012).


