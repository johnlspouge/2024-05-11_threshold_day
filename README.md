# 2024-05-11_threshold_day/

Contains 3 directories:
1. **Data/**
2. **Executable/**
3. **Output/**

### Structure of Executable/ ###

**Executable/** contains 5 subdirectories:

0. modules/ : reusable subroutines that required by the executables below.
modules/ contains a test file test_make.py.
> python test_make.py

tests the installation and functionalities of the Python files. If all is well, the program exits without complaint.

The remaining subdirectories contain executables with specific functions. 

1. Exponent/
2. Extinction_Probability/
3. Quantile/
4. Threshold/

Most of the subdirectories for executables have analogous files, as follows. 

1. run_ui_[executable].py : A program built with the Python argparse package, so '-h' displays program arguments. 
2. run_ui_[executable]_make.py : A Makefile that runs the program.
<br>The Makefile therefore contains an example of the Python command to run the file, along with examples of the input and output filenames relevant to the program.
4. run_ui_[executable].py : A text file, sometimes blank, that received the print statements from the Python program.
<br>The text file often indicates in a long running program where a failure occurred or how long execution took.
5. requirements.txt : contains the version numbers of the Python packages used by the program.

The executables in the subdirectories of **Executable/** have the following functions, also discoverable by running the Python program with the '-h' option.

1. Exponent/run_ui_exponent.py :
<br> Calculates the exponential growth lambda for SEIR model, where E and I are gamma-
distributed.
2. Extinction_Probability/run_ui_extinction_probability.py :
<br> Calculates the extinction probability for (sufficiently simple) single-type Galton-
Watson processes.
3. Quantile/run_quantile_confidence_interval.py : (note the absence of the 'ui_' substring in the filename)
<br> Calculates the quantiles of an integer distribution.
   1. The input is a Dataframe with the arbitrary left-padding columns. After the left-padding columns, columns have headings, strings corresponding to consecutive integers. The columns underneath the headings contain counting numbers or 0. The numbers represent, e.g., the output of a Monte Carlo program that counts the realizations producing a particular integer value, e.g., a particular threshold day (see below).
   2. The output is a Dataframe that first copies the arbitrary left-padding columns. After the left-padding columns, columns have headings corresponding to the quantiles input into the program. Each quantile (e.g., '0.025', '0.05', etc.) corresponds to 3 columns: 'left', 'sample', and 'right'.  The entries in the 'sample' column give the actual quantile in the input. The 'left' and 'right' columns give the endpoints of a bootstrapped confidence interval, to indicate the stability of the sample quantile. 
4. Threshold/run_ui_threshold.py : 
<br> Calculates the threshold day in direct simulation of a superspreading SEIR compartmental model with imperfect ascertainment. It does not use a branching process approximation, which explains the complexity of computation below. The remainder of the README.md explains the simulation for the threshold day as a special case of a very general epidemiological simulation program.

Time may be either continuous or discrete with a constant time-step. Epidemic compartments are indexed and the code tracks their totals (but not their ancestral relationships) as time progresses. Each possible movement from one compartment to another is modeled as a 2-tuple, a "transit". Changes at a given time-step are modeled as a dictionary of "transitions", whose keys are transits and whose values are the counts of the individuals making the transit. A compartment may be stable with no transits (e.g., "Recovered" in SEIR models); others, stable except for susceptibility to transits corresponding to infection (e.g., "Susceptible"); others, subject to spontaneous transits called "decays" (e.g., "Exposed"); and yet others, infectious (e.g., "Infectious"), because they cause transits from one compartment to another (e.g., from "Susceptible" to "Exposed"). The program permit stratification of a population by, e.g., age or susceptibility.

The code is specifically designed to handle non-Markovian decays, because they can cause super-spreading in epidemics [1]. A dictionary of "events", whose keys are times and whose values are the corresponding transitions, holds future non-Markovian decays. A heap sorts the dictionary keys, the times, so the next non-Markovian time can be determined, and if desired popped from the heap-dictionary. 

The computation is initialized with the compartments, their totals, their possible decays, and the infectious transits they induce, if any. The heap-dictionary holds the future decay events for all compartments. The time of the earliest decay is noted. Up to that time but before any infection occurs, infection is modeled as a memoryless process, so the current time plus an exponential or geometric variate is compared to the earliest decay. If the infection is not later than the decay, the infection takes place and is pushed onto the heap-dictionary and then immediately popped with any simultaneous decays. Otherwise, the infection is ignored (the memoryless process, truncated) and the heap-dictionary is popped. In discrete time, therefore, infection takes place before any simultaneous decay. After popping the heap-dictionary, the compartment totals and the current time are updated. 

The continuous-time computation simulates a standard probabilistic epidemic model. It can be overwhelmed, however, by making the infected population large enough, because its computational time is proportional to the infected population size over the modeled duration of the epidemic. The discrete-time version discretizes the computation so the counts in the transitions, not the computational time, reflect the population size. Moreover, as the discrete-time step-size approaches 0, the continuous-time computation is its limit. With a particular view to modeling the early epidemic accurately, however, the discrete-time version can match the initial branching process distribution of an epidemic if the transmissibility &beta; is constant. A deterministic continuous recovery time gives a Poisson individual basic reproduction number, the same as a deterministic discretized recovery time with a Poisson distribution of infection at each time-point. A gamma-distributed continuous recovery time gives a negative binomial individual basic reproduction number, the same as a discretized recovery time with a negative binomial duration and a binomial distribution of infection at each time-point.

The code's initial match to branching processes has been verified by comparing the simulation of SIR and SEIR epidemics from a single infectious individual to analytic formulas for extinction probabilities. It has also been verified by exploiting the flexibility of its compartmental model, as follows. Start an "epidemic" with an single Infectious individual who decays to an inert Recovered. Before its decay, the lone Infectious individual infects Susceptibles, which the simulation converts into inert Exposeds. After repeatedly simulating the decay of the lone Infectious individual, the Exposed count was matched to the distribution of the individual basic reproduction number.

The computation requires a restriction that Infectious compartments can only infect stable compartments. The restriction means that the computation pops only the head decay from the heap-dictionary, whereas without the restriction, infection would necessitate popping decays other than the head decay.
