This project is based on the article "Federated Survival Analysis with Discrete-Time Cox Models" by Andreux et al. - 2020. (https://arxiv.org/pdf/2006.08997)

This article proposes using a discrete-time model as a way to bypass the non-separability of the loss function. This method is also interesting because very little data is exchanged between sites and confidentiality seems to be respected.

In the code "Algorithm from paper - Adam", the Adam (Adaptive Moment Estimation) optimization method is used to estimate the parameters needed (alpha and beta). This method allows us to find the paramaters regardless of initial values, but needs a lot of communication between sites.

In the code "Algorithm from paper - NR", the Newton-Raphson method is used for estimation. This method needs a lot less communication between sites, but creates big matrices that are slow to compute.
