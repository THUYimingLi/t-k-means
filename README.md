This is the implementation of our paper [t-k-means: A Robust and Stable k-means Variant](https://arxiv.org/pdf/1907.07442.pdf), accepted by the IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2021. The project is developed based on the MATLAB, created by [Yang Zhang]<seednov@outlook.com> and [Yiming Li](http://liyiming.tech/). In this paper, we propose a novel robust and stable k-means variant, the t-k-means, and its fast version based on the understanding of k-means, GMM, and TMM.



# Citation
If our work is useful for your research, please cite our paper as follows:

```
@inproceedings{li2021t,
  title={t-k-means: A Robust and Stable k-means Variant}
  author={Li, Yiming and Zhang, Yang and Tang, Qingtao and Huang, Weipeng and Jiang, Yong and Xia, Shu-Tao},
  booktitle={ICASSP},
  year={2021}
}
```

# Description of Main Codes
* **main.m**: test algorithms on a specified dataset and generate results.
* **gmmCluster.m**: the implement of GMM algorithm.
* **tmmCluster.m**: the implement of TMM algorithm.
* **kmeansCluster.m**: the implement of k-means algorithm.
* **kmeansppCluster.m**: the implement of k-means++ algorithm.
* **kmedianCluster.m**: the implement of k-median algorithm.
* **kmedoidCluster.m**: the implement of k-medoid algorithm.
* **SigmaAlphaCluster.m**: the implement of t-k-means algorithm.
* **Sigma0Cluster.m**: the implement of fast t-k-means algorithm.
* **Sigma0ppCluster.m**: the implement of fast t-k-means++ algorithm.
* **plot_loss.m**: loss visualization.

# Evaluation
run 
```
main.m
```
with different parameters (i.e., data_version, method_count = 2, exp_count = 10) and method names to generate needed results.


