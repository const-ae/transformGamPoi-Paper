#!/bin/bash

if [[ $# -eq 0 ]] ; then
  for notebook in notebooks/*.Rmd; do
    R -e "rmarkdown::render('$notebook', output_format = 'html_document')"
  done
else
  for notebook in $@; do
    R -e "rmarkdown::render('$notebook', output_format = 'html_document')"
  done
fi

wait
echo "Done rendering"

rm notebooks/*-tikzDictionary

echo "Done"