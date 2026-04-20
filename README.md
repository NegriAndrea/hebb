# hebb
Halo Extreme Bayesian Bootstrap package to estimate the distribution of N heaviest halos masses you can find in a survery at a given redshift.

`hebb` uses the Uchuu simulation halo catalogue (https://skiesanduniverses.org/Simulations/Uchuu) by sampling sub-volumes that match the survey’s effective volume, and recovering the mass distribution of the N largest halos formed at a particular $z$, with an uncertainty estimate. The volume of the survey is computed given the survey's field of view and redshift depth, or manually selected. To validate the output, the code uses a bayesian bootstrap to compute the variance of the inferred parameters.

Optionally, the package can export the list of halos identified in each realization for follow-up analyses, such as tracing halo evolution over time using the Uchuu merger tree.

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

#### Database Setup
In addition to installing the Python package, you must download a reduced version of the Uchuu database. Two versions are available, with a different cut in mass; a light [one](https://uses0-my.sharepoint.com/:u:/g/personal/anegri_us_es/IQCDx5T1XlTnSbM2YmYcFoYMATIk_pI1XKPzp4ycUJ6N12M?e=pDdoDg) (2.3 GB) and a more complete [one](https://uses0-my.sharepoint.com/:u:/g/personal/anegri_us_es/IQBrGP-3e0xHR6pV0YkfCb1hAYS1KLSNnjmDvVb3H6ytsUA?e=fcudjY) (~40 GB). To set the database path you can define the following environment variable in your `~/.bash_profile`
```
# bash
export HEBB_DB_PATH=/path/to/database
```
where you have to change `/path/to/database` to the path of the downloaded database. By default, the code looks first for the most complete catalugue, and if it is not found, it will look for the light one.

NOTE: In order to keep the file size manageble the database contains only log10(M200), positions and merger tree IDs; the halo positions have been binned on a 40 ckpc gridsize and stored as `uint16` integer, which is precise enough for volumes that are usually way larger than 1 cMpc.

#### Test installation
The easiest way to test if everything is correctly installed and set up is to run the command
```
hebb 0. -L 1000 -M 1e14
```

## Quickstart
To estimate the mass distribution of the 5 most massive halo at z=3 (with uncertaintes) for a survey that spans from $z=3$ to $4$, having field-of-view of 100 arcsec²: 
```
hebb 3. --survey 3. 4. 100.  -n 5 -M 1e12 --plot -t
```
where we limited the catalogue galaxies to be larger than $10^{12}~M_\odot$ for a faster database load (see later for usage of `-M`). The option `-t` creates a table with the sampled haloes, their IDs to trace them back and forth in time with the Uchuu merger tree, and a table containing a list of percentiles of the mass distributio for every mass rank. The `--plot` option draws the mass probability distributions.

## Full Usage
The simplest way to use `hebb` is via command line, `hebb -h` returns the user manual:
```
usage: hebb [-h] (--survey z_min z_max fov | -L L) [-n N] [--niter NITER] [--bb BB] [-t] [--plot] [-M M] [--force-light] z_target

Compute the N most massive dark matter haloes that you can find in a given survey with [z_min, z_max] and field-of-view, and their mass distribution, with a non-overlapping or Monte Carlo (overlapping) sub-boxes sampling of the the Uchuu (2/h cGpc)^3 run. Optionally perform a non-parametric bayesian bootstrap to validate the result.

positional arguments:
  z_target              Redshift of your target

options:
  -h, --help            show this help message and exit
  -n N                  Track the N most massive haloes in each box (1 tracks only the most massive) [default: 1]
  --niter NITER         Do a Monte Carlo selection with NITER overlapping boxes instead of the default non overlapping box search [default not active]
  --bb BB               Number of iterations of bayesian bootstrap to estimate the variance of the results [default: 0]
  -t                    Create a table with the sampled haloes
  --plot                Show a plot of the M200 distribution
  -M M                  OPTIMIZATION: Database mass cut, greatly speed up the database loading and tree building but you can incur into empty boxes [default: None]
  --force-light         DEBUG: force the reading of the light catalogue first

Processing mode (required, mutually exclusive):
  --survey z_min z_max fov
                        Survey min z, max z, FOV in arcmin^2 (used to compute box volume)
  -L L                  Size of the box in cMpc, alternative to the boxsize computation from the FOV and z-depth of the survey
```

## API
TBD



## Citations
```
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
```
