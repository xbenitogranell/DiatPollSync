# DiatPollSync

Repository to reproduce analyses from the paper [**'Human practices are behind the aquatic and terrestrial decoupling to climate change in tropical Andes'**](https://www.sciencedirect.com/science/article/abs/pii/S0048969722012074?via%3Dihub) by Xavier Benito, Blas Benito, Maria I Velez, Jorge Salgado, Tobias Schneider, Liviu Giosan & Majoi Nascimento.

Numerical analyses use recent advances on  [nonlinear estimated rate of change](https://academic.oup.com/icesjms/advance-article/doi/10.1093/icesjms/fsaa056/5835266), GAMs, and [simulated lagged time-series](https://zenodo.org/record/2599859#.YCTzbGhKiUk) to understand time-varying relationships based on assemblage-wide rate of change of aquatic, terrestrial and human datasets using Lake Llaviucu (Ecuador) as a model system.

### Workflow
The Rmd document shows the most relevant code used to reproduce analyses and generate outputs. See R files within the folder "scripts" for a more specific code. Please refer to the full manuscript for full interpretation and data and plots.


### Repository structure
1. **R scripts**: codes for running all the statistical analyses in a logical order
2. **Data**: raw datasets (diatom sediment and training set abundances, lakes' environmental data, and geochemical proxies), and secondary data matrices resulting from primary statistical analyses in <i>R scripts</i>
3. **Figures**: plots resulting from analyses in <i>R scripts</i>

### Graphical abstract
 <img src="outputs/graphical abstract_STOTEN_R1.png" width=600></img>
