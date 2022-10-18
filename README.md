# Calibrating Voronoi cell-based model by SMC-ABC and SNL
![Inference](https://img.shields.io/badge/Inference-Simulation--based%20inference-green)
![Simulation](https://img.shields.io/badge/Simulation-Agent--based%20models-blue)


This repo contains a Python and C++ implementation for following papers:
- [Calibration of a Voronoi cell-based model for tumour growth using approximate Bayesian computation](https://www.biorxiv.org/content/10.1101/2022.09.13.507714v3.abstract)
- Being Bayesian in the 2020s: opportunities and challenges in the practice of modern applied Bayesian statistics

by Xiaoyu Wang, [Adrianne L. Jenner](https://www.adriannejenner.com/), [Robert Salomone](https://robsalomone.com/) and [Christopher Drovandi](https://chrisdrovandi.weebly.com/)

The code to simulate the VCBM originates from [Examining the efficacy of localised gemcitabine therapy for the treatment of pancreatic cancer using a hybrid agent-based model](https://www.biorxiv.org/content/10.1101/2022.04.18.488716v1.abstract)
***
Agent-based models are a class of models that can describe complicated phenomena by analysing the interactions between each agent. In cellular dynamics, treating each cell as an agent is a popular strategy. Prior research has sometimes ignored the geometries of agents and only considered their interactions. In the Voronoi cell-based model, the agent is modelled with Voronoi tesselation so that its shape is as accurate as reality.

In this work, we use sequential Monte Carlo - approximate Bayeisan computation (SMC-ABC) and sequential neural likelihood method to quantify the uncertainty of the model paramters. 


## What does this code do?






## Tips



## Reference
If you find the code useful for your research, please consider citing
- For VCBM
```bib
    @article{jenner2022examining,
    title={Examining the efficacy of localised gemcitabine therapy for the treatment of pancreatic cancer using a hybrid agent-based model},
    author={Jenner, Adrianne and Kelly, Wayne and Dallaston, Michael and Araujo, Robyn and Parfitt, Isobelle and Steinitz, Dominic and Pooladvand, Pantea and Kim, Peter S and Wade, Samantha J and Vine, Kara L},
    journal={bioRxiv},
    year={2022},
    publisher={Cold Spring Harbor Laboratory}
```

This work is an extension of 
```bib
  @article{wang2022calibration,
  title={Calibration of a Voronoi cell-based model for tumour growth using approximate Bayesian computation},
  author={Wang, Xiaoyu and Jenner, Adrianne L and Salomone, Robert and Drovandi, Chris},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
```

The SNL is powered by 
```bib
  @article{tejero2020sbi,
  title={SBI--A toolkit for simulation-based inference},
  author={Tejero-Cantero, Alvaro and Boelts, Jan and Deistler, Michael and Lueckmann, Jan-Matthis and Durkan, Conor and Goncalves, Pedro J and Greenberg, David S and Macke, Jakob H},
  journal={arXiv preprint arXiv:2007.09114},
  year={2020}
}
```
