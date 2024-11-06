#!/bin/bash

demos=$@
for demo in $demos
do
    (
        echo "Creating Notebook $demo"
        cd $demo
        cat ../colab_header.py run.py plot.py | sed -e "s/DEMO/$demo/g" | jupytext --to notebook --execute -o $demo.ipynb 
    )
done

