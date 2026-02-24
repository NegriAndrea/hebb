# hebb
Halo Extreme Block Bootstrap package to estimate the distribution of heaviest halo mass you can find in a survery at a given redshift.

`hebb` uses as a basedata the Uchuu simulation halo catalogue (https://skiesanduniverses.org/Simulations/Uchuu) to perform a block bootstrap by shooting N boxes of a volume equal of the estimated volume of a survey, and recovering the largest halo formed at a particular z, with an uncertainty estimate. The volume of the survey is computed given the survey's field of view and redshift depth, or manually selected. Optionally, the code can dump the list of halo found in the search, which can be used to perform a trace back in time with the Uchuu merger tree.

Please cite the following papers if you use `hebb` in your work: Negri & Belli (2026), [Ishiyama et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.506.4210I) (Uchuu Data Release 1).

## Installation
There are 2 ways to install `hebb`: via `PyPI`
```
python3 -m pip install hebb
```
Or directly from the repository `https://github.com/NegriAndrea/hebb/` and install locally
```
git clone https://github.com/NegriAndrea/hebb
cd hebb
python3 -m pip install -e .
```

### Database Setup
In addition to installing the Python package, you must download a reduced version of the Uchuu database (~40 GB) from [here](https://uses0-my.sharepoint.com/:u:/g/personal/anegri_us_es/IQBrGP-3e0xHR6pV0YkfCb1hAYS1KLSNnjmDvVb3H6ytsUA?e=fcudjY). To set the database path you can define the following environment variable in your `~/.bash_profile`
```
# bash
export HEBB_DB_PATH=/path/to/database
```
where you have to change `/path/to/database` to the path of the downloaded database.

NOTE: In order to keep the file size manageble the database contains only log10(M200), positions and merger tree IDs; the halo positions have been binned on a 40 ckpc gridsize and stored as `uint16` integer, which is precise enough for volumes that are usually way larger than 1 cMpc.

## Usage
The simplest way to use `hebb` is via command line, `hebb -h` returns the user manual:





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
