# OCBianchiSymbols
An implementation of (overconvergent) Bianchi modular symbols

## Set up
Run the makefile to pre-parse the sage files to produce python files, these can then be included as usual.

## Use
### Classical Modular Symbols
For classical symbols, it is enough to instantiate an instance of `ModSymbSpace`, and then either manually create a modular symbol (`createModularSymbol`), compute one from a basis (`getBasis`), or load a previously saved modular symbol (`loadModularSymbol`). Modular symbols can be called as functions on matrices, this evaluates the symbol on the divisor corresponding to that matrix. They also support arithmetic (addition and scalar multiplication on the left), and Hecke operators (`hecke`). Modular symbols can be saved to file for future use (`save`), which overwrites the file when called.

#### Notes
Currently only levels that are square free (as ideals) are supported, and modular symbols with trivial character. Only the Hecke operators T_ell and U_ell for ell prime can be computed, other operators should be computed via the usual formulae.


### Overconvergent Modular Symbols
Overconvergent modular symbols work much like classical symbols. Spaces of overconvergent modular symbols are instances of `OCModSymbSpace`, and symbols can be created like classical symbols, with one new additional method, `liftClassical`, which lifts a classical symbol (which must be an eigensymbol for U_p). When loading overconvergent symbols are loaded from file, by default, they are assumed to be the second entry in the file (separated from the first by a line of 10 # symbols), this can be overriden with the `appended` optional argument. Similarly, the default for saving a modular symbol is to append it to the file, with a separator, but this can be changed with the `append` optional argument. Overconvergent modular symbols have an additional method, `deformationDirection`, which computes the "slope" of the infinitesimal weight change it uniquely deforms along (for example, a form deforming along parallel weight returns 1). If the point is not smooth, and the deformation direction is not unique, it returns `False`.
