Mqrcm 1.1
=============

Changes with respect to version 1.0
------------------
* require pch >= 1.4
* bug fixed. Replaced class(obj) == "try-error" with inherits(obj, "try-error")
* Replaced pch:::fun with fun <- getFromNamespace(fun, ns = "pch"), added getFromNamespace to imports
* updated my e-mail address
* updated references with issue and page numbers
* fixed warnings with PDF < 0
* fixed predict (contrasts)
* fixed plf with scalar input
* added the sigma part to the loss: this guarantees that larger models have smaller loss