![HPFolding](screenshots/HPFolding.png)

## Introduction

This is an attempt to implement an exact algorithm that solves the protein folding problem in the 2D HP model. 

The Hydrophobic-Polar model is a simplified model proposed by Ken A. Dill in 1985 in which the proteins are represented by a sequence of two types of amino acids: Hydrophobic (H) and Polar (hydrophilic) (P). 

The 2D HP folding problem consists in finding an "optimal" conformation, in a 2D grid, of a protein represented by a sequence of 'H's and 'P's. 

The optimality of a conformation is defined by its score: the number of hydrophobic (H-H) bonds (number of H's that are neighbors on the grid but not in the chain). A conformation is optimal if there are no other conformation for the same protein with a higher score.

The following images show some examples of optimal conformations with their scores.


![ex1](screenshots/ex1_mini.png)
![ex2](screenshots/ex2_mini.png)
![ex3](screenshots/ex3_mini.png)


The aim of this project is to solve this problem using an exhaustive search algorithm. 

## More on the implemented algorithm

This work is based on the MINWALK algorithm described in [this paper](https://www.ncbi.nlm.nih.gov/pubmed/23899013).

The algorithm is basically a generic backtracking algorithm with 6 pruning criteria (those that are described in the paper).


## Project structure

I used Xcode during development primarily for its integrated debugger. However, the whole project can be compiled using the makefile. 

The xCode project contains 4 Targets: 
- HPFoldingDebug: same as the makefile target `hpfold` but with Xcode's debugger
- testGrid: I use it to manually test . Nothing interesting there.
- hpTests: same as `hptest` (see bellow) but with debugger.
- generateTests: the program used to generate the test files (see bellow the hptest section).

### Makefile targets

Here are described in more details, the targets of the makefile.

#### hptest

The main for this target is in [tests/test.cpp](tests/test.cpp).


When launched, it will load the chains contained in the file specified in `testInputsFilePath`, calculate their optimal scores and compare them with the scores contained in the file specified by `testResultsFilePath` .

The structure of those files are 1 chain per line and one max-score per line respectively.

There are 3 pairs of test files already generated: 

- [test_inputs10.txt](tests/test_inputs10.txt) and [test\_expected\_results10.txt](tests/test_expected_results10.txt) contain respectivelly all the possible chains of length 10 and their maximum scores. They're generated with the `generateTest.cpp` program.

- [test_inputs13.txt](tests/test_inputs13.txt) and [test\_expected\_results13.txt](tests/test_expected_results13.txt) are the same but for the chains of length 13.

- [test_inputsRnd.txt](tests/test_inputsRnd.txt) and [test\_expected\_resultsRnd.txt](tests/test_expected_resultsRnd.txt) are the same but with randomly generated chains.


#### hpfold

The main for this target is in [src/main.cpp](src/main.cpp).


#### hpfoldgui

The main is the same as for hpfold, but the the `GRAPH` parameter is passed to gcc when compiling (see the [makefile](makefile)).

This one is basically the same as the previous one but with a graphical interface. 

The UI is drawn using the GLUT framework.


Note that the performances of this version aren't as good as in the CLI-based version.

## Usage

### Compilation
Depending on which version you want to compile, run one of the following commands in the root of the project.


```sh
make 				# console version
make hpfold 		# console version
make hpfoldgui 		# graphical version
make tests 			# tests
```


The resulting executable(s) should appear in the `dist` folder.

### Fold a chain

Just run `./hpfold` or `./hpfoldgui` from the dist directory (depending on the version you want) without any param, you will be prompted to enter the chain you want to fold using only `H` and `P` characters.

Example: 

```sh
$ > ./dist/hpfold
Chain: 	HPPHHPHHHPHHPH
HPPHHPHHHPHHPH		6
```

The max score for this chain is 6.


### Fold several chains

You can pass the path to a file containing a list of chains (one per line) as first argument to hpfold. 

```sh
$ > ./dist/hpfold chains.txt
HHHHHHHHHH	4
PHHHHHHHHH	4
HPHHHHHHHH	4
PPHHHHHHHH	3
HHPHHHHHHH	4
```

If you specify a second file in the command line, it will be used to store the results.


### Run the tests

```sh
$ > ./dist/hptest
```

As explained earlier, you have to specify the input files in [test.cpp](tests/test.cpp).

### Run benchmark

The script [benchmark/run.sh](benchmark/run.sh) runs `hpfold` for each chain in the specified input file (benchmark.txt by default) and measure the execution time with the `time` command.

This is not the most precise benchmarking method since it will measure the total execution time of the command and not only the time spent by the folding process, but it gives a pretty good idea of the performances of the implementation.



# Results & Benchmarks


Here's the execution times for 10 random chains of length 30: 


| **Chain**                        | **Score** | **Time** |
|----------------------------------|-----------|----------|
| `HPPPHHPHPPHHHPPPHHHPPPHHPHPHHH` | 12        | 10.135s  |
| `PPPHHHPPPHPHPHHHHPPPPPHHHHPHPH` | 10        | 0.586s   |
| `HHHHHHPPPHHHHHHPHPHHHPHHPPHHHH` | 17        | 0.589s   |
| `HHHPHPPPHPHHPHPHHPPHPHHPPHPHPP` | 13        | 0.008s   |
| `HPHHPHPPHPPPPHPPPPHHPHPPHHPPHH` | 11        | 2.108s   |
| `HPPHHPHHHPPPHPPPHHHPPPHHHPPHPP` | 10        | 0.711s   |
| `HPHPHPHHHPPPPHHPHPHHPHPHPHPPPP` | 11        | 0.173s   |
| `PPHPHHPPPHHPHPHHHPPHPHPPPPHHHP` | 10        | 1.241s   |
| `PHPHHHPHHPPPHPHPPHHPPPPPHHPPPP` | 8         | 0.349s   |
| `PPHPHHHPPPHPHHHPPHPPPHPHHPHPHH` | 11        | 19.376s  |


Here the result for the standard HP benchmark: 

| **Chain**                                                          | **Length** | **Score** | **Time** |
|--------------------------------------------------------------------|------------|-----------|----------|
| `HPHPPHHPHPPHPHHPPHPH`                                             | 20         | 9         | 0.003s   |
| `PPHPPHHPPPPHHPPPPHHPPPPHH`                                        | 25         | 8         | 0.005s   |
| `PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP`                             | 36         | 14        | 0.013s   |
| `PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPPHHPPHPPHHHHH`                 | 48         | 23        | 9.980s   |
| `HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH` | 64         | 42        | 2.097s   |



# References

- Ken A. Dill. Theory for the folding and stability of globular proteins. Biochemistry, 1985, 24 (6), pp 1501–1509. [http://pubs.acs.org/doi/abs/10.1021/bi00327a032](http://pubs.acs.org/doi/abs/10.1021/bi00327a032)

- Giaquinta E., Pozzi L. An effective exact algorithm and a new upper bound for the number of contacts in the hydrophobic-polar two-dimensional lattice model. J Comput Biol. 2013 Aug;20(8):593-609. [https://www.ncbi.nlm.nih.gov/pubmed/23899013](https://www.ncbi.nlm.nih.gov/pubmed/23899013)



