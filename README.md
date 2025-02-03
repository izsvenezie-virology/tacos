# Tacos
A simple python tool to create plots from coverage data.
Plots can be interactive or printed in a PDF file.

## Installation
To install tacos you need to register to [IZSVenezie's GitLab](https://gitlab.izsvenezie.it).
Open your console and enter:
```
git clone https://gitlab.izsvenezie.it/EdoardoGiussani/tacos
cd tacos
pip install .
```
Finish! Now you can plot your beautiful coverage in beautiful charts.

## Standard plot
Shows the coverage for each region of the genome.
To plot your coverage data just type:
```
tacos -o output.pdf coveragefile.cov
```

## Incremental plots
Shows the number of positions of the genome with a specific coverage.
Option ```-i``` allows to print this kind of chart:

```
tacos -i -o output_incremental.pdf coveragefile.cov
```

## Outputs
When the ```-o``` option is specified the plot will be saved in a file in PDF format. If this option is omitted the plots are displayed in interactive mode.
