# extended P-GPLVM
Forked from LMT (latent manifold tuning model)

The old name of the model is Poisson-Gaussian process latent variable model (P-GPLVM). 

**Description:** Estimate latent variables and tuning curves from spike train using P-GPLVM.

Original P-GPLVM is as shown in
Gaussian process based nonlinear latent structure discovery in multivariate spike train data
[Wu *et al* 2017](http://papers.nips.cc/paper/6941-gaussian-process-based-nonlinear-latent-structure-discovery-in-multivariate-spike-train-data).

Extended P-GPLVM now has functions to infer latent variables of new test data in the learned latent space (same neurons), as described in [Luo et al, 2024](https://rnel.rice.edu/pubs/Luo%20et%20al.%20-%202024%20-%20Extended%20Poisson%20Gaussian-Process%20Latent%20Variable%20.pdf)

Usage
=====

* Launch matlab and cd into the directory containing the code
 (e.g. `cd code/poisson-gplvm/`).

* Examine the demo scripts for annotated example analyses of simulated
datasets: 
	*  `demo1_1DGP.m` - Tutorial script illustrating P-GPLVM for 1-dimensional latent variable with tuning curves generated from 1D Gaussian Process.
	*  `demo2_1DBump.m` - Tutorial script illustrating P-GPLVM for 1-dimensional latent variable with tuning curves generated from 1D Gaussian bumps.
	* `demo3_2DBump.m` - Tutorial script illustrating P-GPLVM for 2-dimensional latent variable with tuning curves generated from 2D Gaussian bumps.

Update
=====
`demo_*_ref.m` are the new demos with a `reference` implementation. The new implementation works for multi-trial data with aligned time points.


