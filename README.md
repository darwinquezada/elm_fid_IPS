<br />
<p align="center"> 
  <h3 align="center">Auto-Encoder Extreme Learning Machine for Fingerprint-Based Positioning: A Good Weight Initialization is Decisive</h3>
</p>

```
                          _____                     
 _________        .------|__o__|------.              
:______.-':      :  .--------------.  :             
| ______  |      | :                : |             
|:______o:|      | |  A-WEAR:       | |             
|:______o:|      | |                | |             
|:______o:|      | |  Loading...    | |             
|         |      | |                | |             
|:_____:  |      | |                | |             
|    ==   |      | :                : |             
|       @ |      :  '--------------'  :             
|       o |      :'---...______...---'              
|       o |-._.-i___/' ._______.  \._              
|'-.____o_|   '-.   '-...______...-'  `-._          
:_________:      `.____________________   `-.___.-. 
                 .'.eeeeeeeeeeeeeeeeee.'.      :___:
               .'.eeeeeeeeeeeeeeeeeeeeee.'.         
              :____________________________:

```


<!-- ABOUT THE PROJECT -->
## Abstract

Indoor positioning based on machine learning models has attracted widespread interest in the last years, given its high performance, usability and potential for data compression. While  Autoencoder Extreme Learning Machine (AE-ELM) has been used to extract meaningful information from the datasets and represent it in a low dimension, k-Nearest Neighbour (k-NN) is still widely used for position estimation. This paper introduces a novel method to initialize the input weights in AE-ELM, namely Factorised Input Data (FID), which is based on the normalized form of the orthogonal component of the input data. We provide a comparative analysis with several traditional ways to initialize the input weights in AE-ELM, showing that FID provides a significantly better reconstruction error. Finally, we have performed an assessment with 13 indoor positioning datasets. The dimensionality of the datasets was reduced more than 11 times on average, while the positioning error suffered a small increment of only 15% (average) in comparison to the baseline.

Authors: Darwin Quezada-Gaibor, Joaquín Torres-Sospedra, Jari Nurmi, Yevgeni Koucheryavy, Joaquín Huerta

### Built With

This framework has been developed using:
* [Matlab](https://www.mathworks.com/products/matlab.html)
* [Python](https://www.python.org/)

## Libraries (python)
* pandas, numpy, seaborn, matplotlib, sklearn, tensorflow


## Datasets 
The datasets can be downloaded either from authors' repository (see README file in original_datasets folder) or from the following repository:

      "Joaquín Torres-Sospedra, Darwin Quezada-Gaibor, Germán Mendoza-Silva,
      Jari Nurmi, Yevgeny Koucheryavy, & Joaquín Huerta. (2020). Supplementary
      Materials for 'New Cluster Selection and Fine-grained Search for k-Means
      Clustering and Wi-Fi Fingerprinting' (1.0).
      Zenodo. https://doi.org/10.5281/zenodo.3751042"


** Copy the datasets (.mat) into **databases** folder.

## Usage

  * To run the experiments please use the following file.

```sh
  Run_EML_FID.m
```

<!-- LICENSE -->
## License

CC By 4.0


<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

The authors gratefully acknowledge funding from the European Union’s Horizon 2020 Research and Innovation programme under the Marie Sk\l{}odowska Curie grant agreement No. $813278$, A-WEAR.
