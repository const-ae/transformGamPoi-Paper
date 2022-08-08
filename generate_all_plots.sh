#!/bin/bash

if [[ $# -eq 0 ]] ; then
  for notebook in notebooks/*.Rmd; do
    if [[ "$notebook" != "notebooks/plot_concept_overview_figure.Rmd" ]]; then
      R -e "rmarkdown::render('$notebook', output_format = 'html_document')"
    fi
  done
  # Make sure that this is run last, so that the necessary is generated
  # before
  R -e "rmarkdown::render('notebooks/plot_concept_overview_figure.Rmd', output_format = 'html_document')"
else
  for notebook in $@; do
    R -e "rmarkdown::render('$notebook', output_format = 'html_document')"
  done
fi

wait
echo "Done rendering"

rm notebooks/*-tikzDictionary

echo "Done"