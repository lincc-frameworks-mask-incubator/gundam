# Gundam 

A package to count pairs at lightspeed and estimate 2-point correlation functions 
of large galaxy samples.

## Getting Started

There is a very decent introduction so you can start using Gundam within 5 minutes
and of course the full documentation. Please read-the-docs [here](https://readthedocs.org/projects/gundam/)

### Prerequisites

You will need to have these:

* [Python 3.5 or later](http://www.python.org/)
* [GCC Compiler (C, Fortran & OpenMP support)](https://gcc.gnu.org/)
* [munch](https://pypi.python.org/pypi/munch)
* [pymorton](https://github.com/trevorprater/pymorton/) (Optional, only needed
if experimenting with different orderings)


### Installing

To install Gundam, you have two choices: (1) build from scratch, or (2) use pip. 
I recommend method (1), since it will allow easy access to modify or extend the 
Fortran counting routines. In any case, make sure to fulfil the required 
dependencies. Option 2 using pip is not yet functional.

You just need to clone the Gundam repository and type make

```
git clone https://github.com/el_samotracio/gundam.git
cd gundam
make
```

By default this will compile and build the library in-place. Feel free to modify 
the Makefile to suit your needs. After compilation, you can optionally install 
the library in your default global-site packages directory

```
python setup.py install
```

or in your default user packages directory

```
python setup.py install --user
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Emilio Donoso** - *Initial work* - [ICATE-Conicet](mailto:emiliodon@gmail.com)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* TODO

