# monoblock_parametric

Scripts for the article [Delaporte-Mathurin, R., Hodille, E., Mougenot, J. et al. Parametric study of hydrogenic inventory in the ITER divertor based on machine learning. Sci Rep 10, 17798 (2020).](https://www.nature.com/articles/s41598-020-74844-w)

1. Run a FEniCS docker container

```
docker run -ti -v $(pwd):/home/fenics/shared --name fenics quay.io/fenicsproject/stable:latest
```
2. Install FESTIM 0.7.1

```
pip install git+https://github.com/RemDelaporteMathurin/FESTIM@0.7.1
```

