# BiSCA: BiSpectral EEG Component Analysis

**Paper:** [The influence of nonlinear resonance on human cortical oscillations](https://doi.org/10.1101/2025.06.27.661950)

[![DOI](https://zenodo.org/badge/1055938308.svg)](https://doi.org/10.5281/zenodo.19643376)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a%2B-orange)](https://www.mathworks.com)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2025.06.27.661950-b31b1b)](https://doi.org/10.1101/2025.06.27.661950)

Official MATLAB implementation of the BiSpectral EEG Component Analysis (BiSCA) framework. This repository contains the code to reproduce the analysis and figures presented in our paper.

## About The Project

A longstanding debate in neuroscience concerns whether macroscale brain signals can be described as purely linear Gaussian processes or they harbor the more complex statistics of nonlinear dynamics. This project introduces **BiSCA**, a framework unifying power spectral and bispectral analysis to test for nonlinearity by identifying inter-frequency harmonic relationships and disambiguating them from superimposed components.

Our key findings show that:
*   The brain's broadband, aperiodic background behaves as a linear, Gaussian process.
*   Narrowband Rho oscillatory components, including the Alpha and Mu rhythms, are the primary source of cortical nonlinearity, exhibiting significant quadratic cross-frequency coupling.
*   There is a clear dissociation between signal power and nonlinearity: while the occipital Alpha rhythm dominates the power spectrum, the strongest nonlinear signatures arise from the less dominant parietal Mu rhythm.

These findings suggest that nonlinear resonance is pervasive in cortical signals, primarily expressed through resonant oscillations rather than aperiodic activity.

## Key Features

*   **Joint Spectral and Bispectral Modeling:** A generalized model that fits the power spectrum and bispectrum simultaneously.
*   **Nonlinearity Detection:** Identifies nonlinear resonance phenomena by modeling harmonically related spectral peaks and their bispectral counterparts.
*   **Component Decomposition:** Additively separates the signal into an aperiodic (linear, Gaussian) **Xi (ξ)** process and an oscillatory (nonlinear) **Rho (ρ)** process.
*   **Statistical Testing:** Provides a robust statistical framework for testing the Gaussianity and linearity of each signal component.

## Getting Started

To get a local copy up and running, follow these simple steps.

### Prerequisites

This code is written in MATLAB. The following toolboxes are required to run the analysis:
*   Signal Processing Toolbox
*   Statistics and Machine Learning Toolbox
*   Optimization Toolbox (required for nonlinear fitting algorithms)

### Installation

1.  Clone the repository:
    ```sh
    git clone https://github.com/rigelfalcon/BiSCA.git
    ```
2.  Open MATLAB, navigate to the cloned directory, and run the `setup.m` script to add all necessary subdirectories to your MATLAB path:
    ```matlab
    setup
    ```

### Usage

To get started, run the `example/BiSCA/demo_BiSCA.m` script located in the `example/BiSCA` directory. This script demonstrates how to apply the BiSCA model to sample EEG data and reproduce key figures from the paper.

```matlab
demo_BiSCA
```

## Reproduce Paper Figures

Scripts in `example/BiSCA/` reproduce the main figures (Fig. 1–4).

### Quick Start

1.  [Download the pre-computed data (~516 MB)](https://1drv.ms/u/c/95d61c7c690aa89f/IQCLAUycyg1AR4-MKj9RIbQ4AQDLsT4Oiah42VzLCf_aI1Y) and extract it under `data/`:
    ```
    data/
      iEEG/           (demo channels)
      sEEG/           (metadata, anatomy)
      EEG/            (head models, scouts)
      result/         (pre-computed BiSCA results)
    ```
2.  Run any figure script:
    ```matlab
    setup                                                          % add paths
    run('example/BiSCA/fig1_bisca_model_fit.m')                    % Fig. 1c-d
    run('example/BiSCA/fig1_spa_decomposition.m')                  % Fig. 1b
    run('example/BiSCA/fig2_nnar_simulation.m')                    % Fig. 2
    run('example/BiSCA/fig3_ieeg_bifrequency_test.m')              % Fig. 3 (iEEG)
    run('example/BiSCA/fig4_component_correlation.m')              % Fig. 4
    ```

### Additional Requirements

*   **GPU (optional but recommended):** Fig. 2 scripts use `gpuArray` for NAR model training; an NVIDIA GPU with the Parallel Computing Toolbox is needed for these.
*   **Figures output** to `figure/` (gitignored).

### Script Index

| Script | Figure |
|--------|--------|
| `fig1_bisca_model_fit.m` | Fig. 1c-d (spectrum + bicoherence fit) |
| `fig1_spa_decomposition.m` | Fig. 1b (SPA decomposition) |
| `fig2_nnar_simulation.m` | Fig. 2g-k (asymmetric NNAR simulation) |
| `fig2_nnar_simulation_sym.m` | Fig. 2b-f (symmetric NNAR simulation) |
| `fig3_eeg_test_topography.m` | Fig. 3a-b (EEG test topographies) |
| `fig3_ieeg_bifrequency_test.m` | Fig. 3d (iEEG bifrequency test) |
| `fig3_ieeg_significance_map.m` | Fig. 3c (iEEG significance map) |
| `fig3_ieeg_peak_frequency.m` | Fig. 3e-f (iEEG peak frequency) |
| `fig4_gaussianity_linearity_violin.m` | Fig. 4a-b (Gaussianity/linearity violin) |
| `fig4_component_correlation.m` | Fig. 4c-f (component correlations) |
| `fig4_spatial_distribution.m` | Fig. 4g-l (spatial distribution) |

## Data Availability

The scalp EEG and intracranial EEG (iEEG) data analyzed in this study were obtained from the following publications:
*   **EEG:** Li et al. (2022), "Harmonized-Multinational qEEG Norms (HarMNqEEG)". *NeuroImage*.
*   **iEEG:** Frauscher et al. (2018), "Atlas of the normal intracranial electroencephalogram". *Brain*.

Please refer to the original publications for information on accessing the datasets.

## Citation

If you use this code or our findings in your research, please cite our paper:

Wang, Y., Li, M., Garcia Reyes, R., Bringas-Vega, M. L., Minati, L., Breakspear, M., & Valdes-Sosa, P. A. (2025). The influence of nonlinear resonance on human cortical oscillations. *bioRxiv*. https://doi.org/10.1101/2025.06.27.661950

```bibtex
@article{Wang2025Nonlinear,
  title={The influence of nonlinear resonance on human cortical oscillations},
  author={Wang, Ying and Li, Min and Garc{\'i}a Reyes, Ronaldo and Bringas-Vega, Maria L. and Minati, Ludovico and Breakspear, Michael and Valdes-Sosa, Pedro A.},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.06.27.661950},
  publisher={Cold Spring Harbor Laboratory}
}
```

## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.

## Contact

Pedro A. Valdes-Sosa - pedro.valdes@neuroinformatics-collaboratory.org

Project Link: [https://github.com/rigelfalcon/BiSCA](https://github.com/rigelfalcon/BiSCA)