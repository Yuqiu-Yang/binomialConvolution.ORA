# Inference for Error-Prone Count Data: Estimation under a Binomial Convolution Framework

![Logo](./assets/logo.jpg)

## Introduction 
Measurement error in count data is common but underexplored in the literature, particularly in contexts where observed scores are bounded and arise from discrete scoring processes. Motivated by applications in oral reading fluency assessment, we propose a binomial convolution framework that extends binary misclassification models to settings where only the aggregate number of correct responses is observed, and errors may involve both overcounting and undercounting the number of events. The model accommodates distinct true positive and true negative accuracy rates and preserves the bounded nature of the data. 

Assuming the availability of both contaminated and error-free scores on a subset of items, we develop and compare three estimation strategies: maximum likelihood estimation (MLE), linear regression, and generalized method of moments (GMM). Extensive simulations show that MLE is most accurate when the model is correctly specified but is computationally intensive and less robust to misspecification. Regression is simple and stable but less precise, while GMM offers a compromise in model dependence, though it is sensitive to outliers. 

In practice, this framework supports improved inference in unsupervised settings where contaminated scores serve as inputs to downstream analyses. By quantifying accuracy rates, the model enables score corrections even when no specific outcome is yet defined. We demonstrate its utility using real oral reading fluency data, comparing human and AI-generated scores. Findings highlight the practical implications of estimator choice and underscore the importance of explicitly modeling asymmetric measurement error in count data. 

## Before you start

To carry out MLE for a Binomial Convolution distribution, you will first need to install our `binomialConvolution` R package. You can find the details [here](https://github.com/Yuqiu-Yang/binomialConvolution).


## `code` folder

The `code` folder contains the code snippets to reproduce the results in our paper. The data we generated are accessible on [Zenodo](https://doi.org/10.5281/zenodo.15750744).

In this folder, you will find several `.R` files. They are repeated called to either simulate passage data or to carry out estimation using different methods.

### `code/passage_simulation` subfolder

In this subfolder, there are three `.R` files used to simulate passages under different settings.

### `code/simulation_est` subfolder

For the three scenarios, we created three subfolders `rmse`, `rmse_misspec`, and `se` to contain the code needed for estimating the coefficients in our model with `gmm`, `mle`, and `reg`.

### `code/real_data_est` subfolder

Similar to the `code/simulation_est`, this one contains the code for estimating the coefficients in our model with `gmm`, `mle`, and `reg` with the real data.

### `code/analysis` subfolder

These scripts will take the estimations generated from the previous two steps and compute various metrics. The filenames should be self-explanatory enough. 


## Citation :book:

```
@misc{yang2025inferenceerrorpronecountdata,
      title={Inference for Error-Prone Count Data: Estimation under a Binomial Convolution Framework}, 
      author={Yuqiu Yang and Christina Vu and Cornelis J. Potgieter and Xinlei Wang and Akihito Kamata},
      year={2025},
      eprint={2506.20596},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2506.20596}, 
}
```


