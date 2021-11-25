BASC: Bootstrap Analysis of Stable Clusters
===========================================
-   [Dataset of Cognitive scores](#dataset-of-cognitive-scores)
-   [Cluster selection](#cluster-selection)
-   [Similarity Matrix](#similarity-matrix)
-   [Similarity Matrix of the bootstraped data](#similarity-matrix-of-the-bootstraped-data)
-   [Stablility Matrix: 10,000 permutations](#stablility-matrix-10000-permutations)
-   [Hierarchical agglomerative clustering apply to the Stability matrix](#hierarchical-agglomerative-clustering-apply-to-the-stability-matrix)
-   [Stable Cluster Matrix](#stable-cluster-matrix)

Dataset of Cognitive scores
---------------------------

Each row represents a subject and the columns the different scales. The dendrograms is based on a hierarchical clusteritation with euclidean distances.
Patient 27872 was excluded from this cluster analysis due to his high overal IQ score.  
&gt; *y*=data set input  
&gt; *P*=partition  
&gt; *Φ*= cluster operation  
  
<img src=
"https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AY+%5Cxrightarrow%7Bf%7D+y+%5Cxrightarrow%7B%5CPhi%7D%5CPhi+y%0A%5Cend%7Balign%2A%7D" 
alt="\begin{align*}
Y \xrightarrow{f} y \xrightarrow{\Phi}\Phi y
\end{align*}">

  
 ![](BASC_files/figure-markdown_github/unnamed-chunk-3-1.png)

Cluster selection
-----------------

based on 26 indexes of clustering, according to the majority rule, the best number of clusters will be selected. NbClusters: **"An `R` Package for Determining the Relevant Number of Clusters in a Data Set"**.

![](BASC_files/figure-markdown_github/unnamed-chunk-4-1.png)

Similarity Matrix
-----------------

Quantifies the stable festures of the stochastic clustering process. Is a matrix where *Φ*<sub>*ij*</sub> = 1 if they are both in the same cluster and 0 in the contrary.  
 
 <img src=
"https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0AS_%7Bij%7D%3D+Pr++%26+%5Cleft%28++%5CPhi_%7Bij%7D%28y%29%3D+1+%5Cmid+Y+%5Cxrightarrow%7Bf%7D+y+%5Cright%29%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
S_{ij}= Pr  & \left(  \Phi_{ij}(y)= 1 \mid Y \xrightarrow{f} y \right)
\end{align*}
">

Stability matrix of the original data, non-bootstrap: ![](BASC_files/figure-markdown_github/unnamed-chunk-5-1.png)

Similarity Matrix of the bootstraped data
-----------------------------------------

Using a monte carlo aproximation we used stratified bootstrap (subjects withs replacement) in order to obtain a *B* independent sample of <img src=
"https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Chat%7Bf%7D%28y%29%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
\hat{f}(y)
\end{align*}
">.  

So the Similarity boot matrix will be the sum of each of the bootstrap similarity matrix.  

<img src=
"https://render.githubusercontent.com/render/math?math=%5CLarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Chat%7BS%7D_%7Bij%7D%3DB%5E%7B-1%7D+%5Csum_%7Bb%3D1%7D%5EB+%5CPhi_%7Bij%7D+%28y%5E%7B%2Ab%7D%29+%3D+%5Chat%7BS%7D_%7Bij%7D%5E%7Bboot%7D+%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
\hat{S}_{ij}=B^{-1} \sum_{b=1}^B \Phi_{ij} (y^{*b}) = \hat{S}_{ij}^{boot} 
\end{align*}
">

Stablility Matrix: 10,000 permutations
--------------------------------------

![](BASC_files/figure-markdown_github/unnamed-chunk-6-1.png)

Hierarchical agglomerative clustering apply to the Stability matrix
-------------------------------------------------------------------

![](BASC_files/figure-markdown_github/unnamed-chunk-7-1.png)

Stable Cluster Matrix
---------------------

![](BASC_files/figure-markdown_github/unnamed-chunk-8-1.png)

Reference
---------------------
For further details in the theory see:
> Bellec, Pierre, et al. "Multi-level bootstrap analysis of stable clusters in resting-state fMRI." Neuroimage 51.3 (2010): 1126-1139.



