# Are shots predictive of soccer results?

## Short description

This case study proposes a multiple outcomes hierarchical model for predicting the results of a soccer match of the English Premier League (EPL) considering the following nested quantities: **number of scores**, **number of shots on target** and **number of total shots**. The model may be seen as a natural extension of the existing ones which only consider the scores (Karlis and Ntzoufras 2003, Baio and Blangiardo 2010) as historical information. Posterior predictive checking, leave-one-out cross validation (LOO) and future predictions are central topics in this presentation.    

## Repository material

- ```egidi_notebook.Rmd```: Rmd file of the notebook.

- ```egidi_notebook.html```: HTML file of the notebook.

- ```egidi_slides.html```: HTML file of the slides from the talk.

- ```PL08-08.csv, PL09-10.csv```, etc.: data files containing scores and shots information about a single EPL season.

- ```shots.stan```: shots Stan model.

- ```double_poisson.stan```: double Poisson Stan model.

- ```Pois_prob.RData```: double Poisson probabilities used for computing the Brier score.

- ```ref.bib```: bibtex list bibliography.


## Computational issues

Due to the dimension of the file, the ```fit.RData``` file containing the Stan fit object is not included in this repository. The interested reader may however reproduce all the code contained in ```egidi_notebook.Rmd``` and save  the Stan fit result (```mod_shots_stan``` in the document) in a new file called   ```fit.RData```.



## Main author contact

Leonardo Egidi

Postdoctoral researcher

Dipartimento di Scienze Economiche, Aziendali, Matematiche e Statistiche

Università degli Studi di Trieste, Trieste, Italy

email: legidi@units.it

