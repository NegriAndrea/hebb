# hebb
Halo Extreme Block Bootstrap package to estimate the largest halo mass given the field of view and redshift depth of a survey. `hebb` uses as a basedata the Uchuu simulation halo catalogue (https://skiesanduniverses.org/Simulations/Uchuu) to perform a block bootstrap by shooting N boxes of a volume equal of the estimated volume of a survey, and recovering the largest halo formed at a particular z, with an uncertainty estimate. Optionally, the code can dump the list of halo found in the search, which can be used to perform a trace back in time with the Uchuu merger tree.
Please cite the following papers if you use `hebb` in your work: Negri & Belli (2026), Ishiyama et al. (2021) (Uchuu Data Release 1).


@ARTICLE{2021MNRAS.506.4210I,
       author = {{Ishiyama}, Tomoaki and {Prada}, Francisco and {Klypin}, Anatoly A. and {Sinha}, Manodeep and {Metcalf}, R. Benton and {Jullo}, Eric and {Altieri}, Bruno and {Cora}, Sof{\'\i}a A. and {Croton}, Darren and {de la Torre}, Sylvain and {Mill{\'a}n-Calero}, David E. and {Oogi}, Taira and {Ruedas}, Jos{\'e} and {Vega-Mart{\'\i}nez}, Cristian A.},
        title = "{The Uchuu simulations: Data Release 1 and dark matter halo concentrations}",
      journal = {\mnras},
         year = 2021,
        month = sep,
       volume = {506},
       number = {3},
        pages = {4210-4231},
          doi = {10.1093/mnras/stab1755},
archivePrefix = {arXiv},
       eprint = {2007.14720},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021MNRAS.506.4210I},
}
