# 2024-05-11_threshold_day/

Contains 3 directories:
1. **Data/**
2. **Executable/**
3. **Output/**

**Executable/** contains 5 subdirectories:

0. modules/ : reusable subroutines that required by the executables below.
modules/ contains a test file test_make.py.
> python test_make.py

tests the installation and functionalities of the Python files. If all is well, the program exits without complaint.

The remaining subdirectories contain executables with specific functions, each containing analogous files. 

1. Exponent/
2. Extinction_Probability/
3. Quantile/
4. Threshold/

The **Data/** and **Output/** directories contain these 4 subdirectories, too, for input and output from the corresponding executables. Each subdirectory for executables generally contains 4 files, as follows. 

1. run_ui_[executable].py : A program built with the Python argparse package, so '-h' displays program arguments.
<br>'-h' displays the purpose of the program (summarized below) and its arguments. 
2. run_ui_[executable]_make.py : A Makefile that runs the program.
<br>The Makefile therefore contains an example of the Python command to run the file, along with the input and output filenames relevant to the program, which need not be in the corresponding **Data/** and **Output/** directories. Each Makefile created an example of program usage, so examining it will yield an example of how to run the program, its arguments, and the input and output filenames for the example.
3. run_ui_[executable].log : A plain text file of system output, sometimes blank. The file received the print statements from the Python program.
<br>In a long running program, the text file often indicates where a failure occurred or how long execution took.
4. requirements.txt : contains the version numbers of the Python packages used by the program.

The subdirectories of **Executable/** contain the following executables with the following purposes, also discoverable by running the Python program with the '-h' option.

1. Exponent/run_ui_exponent.py :
<br> Calculates the exponential growth lambda for an SEIR model where E and I are gamma-
distributed.
2. Extinction_Probability/run_ui_extinction_probability.py :
<br> Calculates the extinction probability for (sufficiently simple) single-type Galton-
Watson processes.
3. Quantile/run_quantile_confidence_interval.py : (note the absence of the 'ui_' substring in the filename)
<br> Calculates the quantiles of an integer distribution.
   1. The input is a Dataframe with the arbitrary columns padding the left. After the left-padding columns, columns have headings, strings corresponding to consecutive integers. The columns underneath the headings contain counting numbers or 0. The numbers represent, e.g., the output of a Monte Carlo program that counts the realizations producing a particular integer value, e.g., a particular threshold day (see below).
   2. The output is a Dataframe. The output copies the arbitrary left-padding columns. After the left-padding columns, columns have headings corresponding to the quantiles input into the program. Each quantile (e.g., '0.025', '0.05', etc.) corresponds to 3 columns: 'left', 'sample', and 'right'.  The entries in the 'sample' column give the actual quantile from the input. The 'left' and 'right' columns give the endpoints of a bootstrapped confidence interval (of user-entered size), to indicate the variability of the sample quantile. 
4. Threshold/run_ui_threshold.py : 
<br> Simulates a general (non-Markovian) superspreading epidemic model, under two restrictions: (1) the force of infection is Markovian (i.e., equivalent to a Poisson process); and (2) only stable compartments become infected. Here, a compartment is "stable" if an individual in it does not change unless infected. In an SIR model, e.g., an individual in the "Susceptible" compartment remains susceptible unless infected. The compartments, age-groups, etc., are otherwise arbitrary.

**A Detailed Explanation of Threshold/run_ui_threshold.py**

A non-Markovian simulation can be computationally expensive. A superspreading model might, e.g., use a general gamma-distributed infectious period. If too many individuals are infected (e.g., greater than 200), an exact non-Markovian simulation becomes computationally burdensome, because each individual's time of recovery from infection must be stored. Two options can reduce computation: (1) an option for approximating with discrete time-steps (not implemented in the present programs, but in principle available with the call set_delta(delta) in the file Executable/modules/jls_epidemic_infect.py); and (2) an option '-w' that creates watchers that terminate the computation when it meets a condition, e.g., if the number of daily infections exceeds y.

The **Threshold/** directories contain data, an executable, and output for calculating a threshold day where 30 or more new cases occur daily. It simulates a superspreading SEIR compartmental model with imperfect ascertainment exactly and uses a watcher to terminate the simulation. It does not use a branching process approximation. 

The actual program can simulate general superspreading epidemics, so its mechanics are explained in the general context.

Time may be either continuous or discrete with a constant time-step. Epidemic compartments are indexed and the code tracks their totals, but not their ancestral relationships, as time progresses. A _transit_ is a 2-tuple of compartments, modeling possible movement from one compartment to another. A _transition_ is a dictionary element whose key is a transit and whose value is a count of individuals making the transit. The corresponding dictionary holds all the transitions that occur at a given time-step. A particular compartment may be _stable_, with no transits in the absence of infection (e.g., "Susceptible" or "Recovered" in SEIR models); or _unstable_ and subject to spontaneous transits called _decays_ (e.g., "Exposed" or "Infectious" in SEIR models). Some compartments are infectious (e.g., "Infectious" in SEIR, but not "Exposed"), because they cause transits from uninfected to infected compartments (e.g., from "Susceptible" to "Exposed"). The program permits stratification of a population by, e.g., age or susceptibility. It also permits random forks to clinical or subclinical infection. 

The code is specifically designed to handle non-Markovian decays, to model super-spreading in epidemics [1]. An _event_ is a dictionary element of whose key is an epidemic time and whose value is a dictionary of decays that occur at the epidemic time. A heap sorts the dictionary keys, the times, so the time of the next decay can be determined, and if desired popped from the heap-dictionary. The algorithm for epidemic simulation ensures that the heap-dictionary of decays holds all future non-Markovian transitions. 

The algorithm is initialized with 3 separate classes of input files: (1) the compartments and their totals; (2) all possible decays; and (3) every possible infectious transit and the infectious compartment that can cause it. The heap-dictionary of decays is initialized with the future decays of all unstable compartments. The earliest time of decay is noted. Up to that time but before any infection occurs, infection is modeled as a memoryless Poisson process. The current time plus the corresponding exponential variate is compared to the earliest time of decay. (In discrete time, the exponential variate is approximated by a geometric variate.) If infection occurs before the decay, it is pushed onto the heap-dictionary, then immediately popped with any simultaneous decays. Otherwise, infection is ignored (the memoryless process is correctly truncated in favor of the earliest decay) and the earliest decay on the heap-dictionary popped. In discrete time, therefore, infection takes place before any simultaneous decay. After popping the heap-dictionary, the compartment totals and the current time are updated. 

The continuous-time computation simulates a standard probabilistic epidemic model. It can be overwhelmed, however, by making the infected population large enough, because its computational time is proportional to the infected population size over the modeled duration of the epidemic. The discrete-time version discretizes the computation so the counts in the transitions, not the computational time, reflect the population size. Moreover, as the discrete-time step-size approaches 0, the continuous-time computation is its limit. With a particular view to modeling the early epidemic accurately, however, the discrete-time version can match the initial branching process distribution of an epidemic if the transmissibility &beta; is constant. A deterministic continuous recovery time gives a Poisson individual basic reproduction number, the same as a deterministic discretized recovery time with a Poisson distribution of infection at each time-point. A gamma-distributed continuous recovery time gives a negative binomial individual basic reproduction number, the same as a discretized recovery time with a negative binomial duration and a binomial distribution of infection at each time-point.

The code's initial match to branching processes has been verified by comparing the simulation of SIR and SEIR epidemics from a single infectious individual to analytic formulas for extinction probabilities. It has also been verified by exploiting the flexibility of its compartmental model, as follows. Start an "epidemic" with an single Infectious individual who decays to an inert Recovered. Before its decay, the lone Infectious individual infects Susceptibles, which the simulation converts into inert Exposeds. After repeatedly simulating the decay of the lone Infectious individual, the Exposed count was matched to the distribution of the individual basic reproduction number.

The computation requires a restriction that Infectious compartments can only infect stable compartments. The restriction means that the computation pops only the head decay from the heap-dictionary, whereas without the restriction, infection would necessitate popping decays other than the head decay.

### Explanation of the Example Data in Threshold/Data/ and Threshold/Output/ ###


