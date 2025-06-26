# Inference for Error-Prone Count Data: Estimation under a Binomial Convolution Framework

![Logo](./assets/logo.jpg)

## Introduction 
Measurement error in count data is common but underexplored in the literature, particularly in contexts where observed scores are bounded and arise from discrete scoring processes. Motivated by applications in oral reading fluency assessment, we propose a binomial convolution framework that extends binary misclassification models to settings where only the aggregate number of correct responses is observed, and errors may involve both overcounting and undercounting the number of events. The model accommodates distinct true positive and true negative accuracy rates and preserves the bounded nature of the data. 

Assuming the availability of both contaminated and error-free scores on a subset of items, we develop and compare three estimation strategies: maximum likelihood estimation (MLE), linear regression, and generalized method of moments (GMM). Extensive simulations show that MLE is most accurate when the model is correctly specified but is computationally intensive and less robust to misspecification. Regression is simple and stable but less precise, while GMM offers a compromise in model dependence, though it is sensitive to outliers. 

In practice, this framework supports improved inference in unsupervised settings where contaminated scores serve as inputs to downstream analyses. By quantifying accuracy rates, the model enables score corrections even when no specific outcome is yet defined. We demonstrate its utility using real oral reading fluency data, comparing human and AI-generated scores. Findings highlight the practical implications of estimator choice and underscore the importance of explicitly modeling asymmetric measurement error in count data. 

## Repo structure 



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


