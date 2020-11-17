# Meta-d_STAN
The Stan codes of Meta-d' Model

- Term Project of Computational modeling course (Spring 2020)
- For the details of the term project, please visit [here](https://drive.google.com/drive/folders/1-pREh7Nv_k1mXIbE79J3DAAJm_To5O0Z?usp=sharing) (Google Drive). You can find a pdf file and other versions of codes.
- The codes implements the meta-d' model ([Maniscalco & Lau 2012](https://doi.org/10.1016/j.concog.2011.09.021)) in a non-hierarchical and hierarchical Bayesian framework using STAN. The hierarchical version also referred HMeta-d' toolbox from [Fleming 2017](https://doi.org/10.1093/nc/nix007).

## Folders and Files
- `/codes`: Stan codes and R codes for running the model. 'indiv' for non-hierarchical Bayesian modeling, 'hierarchical' for hierarchical Bayesian modeling
- `/example_data`: Example data for testing the codes. The processed open data from [Samaha & Postle 2017](http://dx.doi.org/10.1098/rspb.2017.2035)
- `/figure`: Example figures

## Reference

- Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.[Paper](https://doi.org/10.1016/j.concog.2011.09.021) doi:10.1016/j.concog.2011.09.021 [MATLAB code](http://www.columbia.edu/~bsm2105/type2sdt/)

- Fleming, S. M. (2017). HMeta-d: hierarchical Bayesian estimation of metacognitive efficiency from confidence ratings. Neuroscience of Consciousness, 3(1), nix007-. doi:10.1093/nc/nix007 [Paper](https://doi.org/10.1093/nc/nix007) [HMeta-d toolbox](https://github.com/metacoglab/HMeta-d)

- Samaha, J., & Postle, B. R. (2017). Correlated individual differences suggest a common mechanism underlying metacognition in visual perception and visual short-term memory. Proceedings of the Royal Society B, 284(1867), 20172035. doi:10.1098/rspb.2017.2035 [Paper](http://dx.doi.org/10.1098/rspb.2017.2035)
