# Optical coherence refraction tomography

*Reproduction de la methode OCRT pour le projet de 4eme annee de baccalaureat Polytechnique Montreal.*
* **Papier original**: https://www.nature.com/articles/s41566-019-0508-1
* **Nature Supplementary Materials**: [PDF](https://static-content.springer.com/esm/art%3A10.1038%2Fs41566-019-0508-1/MediaObjects/41566_2019_508_MOESM1_ESM.pdf)
* **Code original**: https://github.com/kevinczhou/optical-coherence-refraction-tomography

## Description algorithme
![alt text](data/description_algorithme.PNG)

## Data
Les images sont trop grosses pour etre directement telecharger sur Github (limite > 100Mb)

## Comment debuter
1) Clone this GitHub repository on your computer:
```bash
git clone https://github.com/PaulBautin/optical-coherence-refraction-tomography.git
cd optical-coherence-refraction-tomography
```
- For Windows user, you might need to [install git](https://git-scm.com/downloads) prior to clone the repository.

2) Run the following to create a virtual environment and start the notebook:

```bash
conda env create -f environment.yml # Only do it once in order to create the environment (might take a few minutes)
```

- Make sure that your prompt is currently on the `optical-coherence-refraction-tomography` folder when you call the `environment.yml` file.
- For Windows user, you might need to type these commands in `Anaconda Prompt` if `cmd` does not recognize `conda`.

**Make sure that you have the last version of the files by pulling the repo before every new lab** (`git pull`).

## Execution
Le code a ete execute avec python 3.6. Voici les informations importantes pour faire tourner le script principal *image_correction.py*
```bash
usage: image_correction [-h] -i I [-fig] [-o O] [-l]

Algorithme de correction des distorsions liées aux changements d'indices de réfractions.

optional arguments:
  -h, --help  show this help message and exit

MANDATORY ARGUMENTS:
  -i I        Path to folder that contains input images. Example: "octr_data"

OPTIONAL ARGUMENTS:
  -fig        Generate figures
  -o O        Path where figures will be saved. By default, they will be saved in the current directory.
  -l          manually find borders of capillary
```
