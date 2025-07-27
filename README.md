TA-PDC: Provable Data Contribution with Traceable Anonymous for Group Transactions

This repository implements TA-PDC, a privacy-preserving provable data contribution scheme with traceable anonymity, supporting fair and efficient quantification of user contributions.

Our code is based on the excellent work [wts](https://github.com/sourav1547/wts) and has been modified and extended accordingly.

Runtime Environment: Implemented in Go 1.19 and use the BLS12-381 asymmetric pairing based curve implementation from gnark-crypto, where each F element is 256 bits, each G1 element is 384 bits, and each G2 element is 768 bits. Simulation experiments are run on Ubuntu 18.04 VMware 17.5 on a laptop running Windows 11 with Intel(R) Core(TM) Ultra 5 125H 3.60GHz 16 GB RAM. 

The source code is located in the src directory with the following structure:

  pdc.go: Core implementation of the TA-PDC scheme.

  pdc_test.go: Correctness verification and performance benchmarking using Go’s built-in testing and benchmarking framework.

  oruta.go: Core implementation of the Oruta: Privacy-Preserving Public Auditing for Shared Data in the Cloud scheme.

  oruta_test.go: Correctness verification and performance benchmarking for the Oruta scheme.

  utils.go & utils_test.go: Utility functions and their corresponding tests, providing auxiliary support for the above schemes.
