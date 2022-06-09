# Continental-scale empirical evidence for relationships between fire response strategies and fire frequency
The repository containing R code and data for analysing relationships between fire response strategies and fire frequency in Australia. 

Authors:
Sophie Yang, Mark K. J. Ooi, Daniel S. Falster, Will K. Cornwell


Methods for extracting fire response data on resprouting and seeding from [AusTraits](https://austraits.org/) are available in "code/extracting_fire_response_data.Rmd".

Methods for calculating the mean FRI across the distribution of each species/taxon using Poisson Regression are available in "code/mean_fri_poisson_regression_method.Rmd". Methods for using the survival analysis method closely following Simpson et al. (2021) are in "code/median_fri_survival_analysis_method".

The methods for analysing relationships between fire response strategies and fire frequency are available in "code/fire_response_strategies_and_fire_frequency.Rmd". The methods for analysising relationships between fire response strategies and leaf traits are available in "fire_response_strategies_and_leaf_traits.Rmd".

The processes of cleaning GBIF data, making an Australian native plant lookup and a lookup for GBIF/APC taxonomy were borrowed from the following repository: https://github.com/traitecoevo/aus_gbif_clean.


References:
**Simpson KJ, Jardine EC, Archibald S, Forrestel EJ, Lehmann CER, Thomas GH, Osborne CP. 2021.** Resprouting grasses are associated with less frequent fire than seeders. *New Phytologist* **230**: 832–844.

