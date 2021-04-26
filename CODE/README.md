# Exploring Gene Expresison Levels in Wound Healing  
This is the self-contained project folder. To reproduce our results you need only this folder.  

Full preprocessed dataset is already provided in this folder:  
But here's a link to it in the GitHub repo we made for this project anyway:  
https://github.com/ZovcIfzm/Kimball-ATVG-outlier-analysis/blob/main/results/count_output.txt  

While we've written a makefile for the project, (run "make run") the code itself is just a python file, which you can run with the command:  

``python3 clustering_outliers.py --counts_file count_output.txt --show_plots 1``  

The only other step needed is to download the python packages which can be done with  

``pip install -r requirements.txt``  