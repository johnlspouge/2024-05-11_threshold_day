# 2024-05-11_threshold_day/

Contains 3 directories:
1. **Data/**
2. **Executable/**
3. **Output/**

**Executable/** contains 5 subdirectories:

0. modules/ : reusable subroutines that the executables below require.
modules/ contains a test file test_make.py.
> python test_make.py

tests the installation and functionalities of the Python files in modules/. If all is well, the program exits without complaint.

The remaining subdirectories contain executables with specific functions, each containing analogous files. 

1. Exponent/
2. Extinction_Probability/
3. Quantile/
4. Threshold/

The **Data/** and **Output/** directories contain these 4 subdirectories, too, so they can contain the input and output from the corresponding executables if desired. Each subdirectory for executables generally contains 4 files, as follows. 

1. run_ui_[executable].py : A program built with the Python argparse package, so '-h' displays program arguments.
<br>'-h' displays the purpose of the program (also summarized below) and its arguments. 
2. run_ui_[executable]_make.py : A Makefile that runs the program.
<br>The Makefile therefore contains an example of the Python command to run the file, along with the input and output filenames relevant to the program, which need not be in the corresponding **Data/** and **Output/** directories. Each Makefile provides an example of program usage, so examining it will yield an example of how to run the program, the program arguments, and the input and output filenames for the example.
3. run_ui_[executable].log : A plain text file of system output, sometimes blank. The file received the print statements from the Python program.
<br>In a long running program, the text file often indicates where a failure occurred or how long execution took.
4. requirements.txt : contains the version numbers of the Python packages used by the program.

The subdirectories of **Executable/** contain the following executables with the following purposes, also discoverable by running the Python program with the '-h' option.

1. Exponent/run_ui_exponent.py :
<br> Calculates the exponential growth lambda for an SEIR model where E and I are gamma-distributed.
2. Extinction_Probability/run_ui_extinction_probability.py :
<br> Calculates the extinction probability for (sufficiently simple) single-type Galton-Watson processes.
3. Quantile/run_quantile_confidence_interval.py : (note the absence of the 'ui_' substring in the filename)
<br> Calculates the quantiles of an integer distribution.
   1. The input is a Dataframe with the arbitrary columns padding the left. After the left-padding columns, columns have headings, strings corresponding to consecutive integers. The columns underneath the headings contain counting numbers or 0. The numbers represent, e.g., the output of a Monte Carlo program that counts the realizations producing a particular integer value, e.g., a particular threshold day (see below).
   2. The output is a Dataframe. The output copies the arbitrary left-padding columns. After the left-padding columns, columns have headings corresponding to the quantiles input into the program. Each quantile (e.g., '0.025', '0.05', etc.) corresponds to 3 columns: 'left', 'sample', and 'right'.  The entries in the 'sample' column give the actual quantile from the input. The 'left' and 'right' columns give the endpoints of a bootstrapped confidence interval (of user-entered size), to indicate the variability of the sample quantile. 
4. Threshold/run_ui_threshold.py : 
<br> Simulates a general (non-Markovian) superspreading epidemic model, under two restrictions: (1) the force of infection is Markovian (i.e., equivalent to a Poisson process); and (2) only stable compartments become infected. Here, a compartment is _stable_ if an individual in it does not change unless infected. In an SIR model, e.g., individuals in the "Susceptible" and "Recovered" compartments remain there unless infected. The compartments, age-groups, etc., are otherwise arbitrary.

**A Detailed Explanation of Threshold/run_ui_threshold.py**

A non-Markovian epidemic simulation can be computationally expensive. A superspreading model might, e.g., use a general gamma-distributed infectious period. If too many individuals are infected (e.g., greater than 200), an exact non-Markovian simulation becomes computationally burdensome, because each individual's time of recovery from infection must be stored separately. Two options can reduce computation: (1) an option for approximating with discrete time-steps (not implemented in the present programs, but in principle available with the call set_delta(delta) in the file Executable/modules/jls_epidemic_infect.py); and (2) an option '-w' that creates watchers that terminate the computation when it meets a condition, e.g., if the number of daily infections exceeds y.

The **Threshold/** directories contain data, an executable, and output for calculating a threshold day where 30 or more new cases occur daily. It simulates a superspreading SEIR compartmental model with imperfect ascertainment exactly and uses a watcher to terminate the simulation. It does not use a branching process approximation. 

The actual program can simulate general superspreading epidemics, so its mechanics are explained in the general context.

Time may be either continuous or discrete with a constant time-step. Epidemic compartments are indexed and the code tracks their totals as time progresses, but not their ancestral relationships. A _transit_ is a 2-tuple of compartments, modeling possible movement from one compartment to another. A _transition_ is a dictionary element whose key is a transit and whose value is a count of individuals making the transit. The corresponding dictionary holds all the transitions that occur at a given time-step. A particular compartment may be _stable_, with no transits in the absence of infection (e.g., "Susceptible" or "Recovered" in SEIR models); or _unstable_ and subject to spontaneous transits called _decays_ (e.g., "Exposed" or "Infectious" in SEIR models). Some compartments are infectious (e.g., "Infectious" in SEIR, but not "Exposed"), because they cause transits from uninfected to infected compartments (e.g., from "Susceptible" to "Exposed"). The program permits stratification of a population by, e.g., age or susceptibility. It also permits random forks to a clinical or subclinical infection. 

The code is specifically designed to handle non-Markovian decays, to model super-spreading in epidemics [1]. An _event_ is a dictionary element of whose key is an epidemic time and whose value is a dictionary of decays that occur at the epidemic time. A heap sorts the time-keys so the time of the earliest decay can be determined, and if desired popped from the heap-dictionary. The algorithm for epidemic simulation (next) ensures that the heap-dictionary of decay events holds all future non-Markovian events that are contingent on the present or past. 

The algorithm is initialized with 3 separate classes of input files: (1) the compartments and their totals; (2) all possible decays; and (3) every possible infectious transit and the infectious compartment that can cause it. The initialization does not permit stable infectious compartments. The computation also requires that Infectious compartments only infect individuals from stable compartments. Otherwise, the infection might necessitates popping a future decay already stored in the heap-dictionary, and not necessarily the earliest decay. The heap-dictionary of decay events is initialized with the future decays of all unstable compartments. The current epidemic time set to 0, and the time of the earliest decay is noted. Infection is then modeled as a memoryless Poisson process. The current epidemic time plus the corresponding exponential variate for the first infection is compared to the earliest epidemic time of decay. (In discrete time, the exponential variate is approximated by a geometric variate.) If infection occurs before the decay, the infective event is pushed onto the heap-dictionary, then immediately popped with any simultaneous decay events. The algorithm therefore prioritizes infection over a simultaneous decay. Otherwise, the algorithm ignores the possibility of infection (i.e., the memoryless process is correctly truncated in favor of the earliest decay) and the earliest decay on the heap-dictionary popped. When the algorithm pops the heap-dictionary of decays events, it always updates compartment totals and the current epidemic time. The algorithm repeats the foregoing steps until the heap-dictionary is empty (no unstable individual remains, so no infectious individual remains), or the watchers have terminated the computation (below).

**Input Files for Executable/Threshold/run_ui_threshold.py**

See **Executable/Threshold/**. For user convenience, its program **run_ui_threshold.py** can simulate several epidemic sequentially. The Makefile **run_ui_threshold_make.py** specifies the following inputs for **run_ui_threshold.py**. The inputs are all in **Data/Threshold/**.  

1. **0_exposed_123.csv** specifies the compartments and their totals. 
<br> The input is a DataFrame with a single column 'compartment' and then several columns 'count'. Column 'compartment' lists all the compartments in the epidemic. Blank cells in the DataFrame are equivalent to 0 entries. Each column 'count' lists the initial compartment count for a different simulation. 
2. **1_add_decay_to_** specifies the possible decays (possibly non-Markovian) for the compartments. 
<br> The input is a DataFrame with 5 columns: 'from', 'to', 'mean', 'dispersion', 'prob'. As one might expect, they specify a decay transit by where it starts and ends, and the mean and dispersion parameters controlling the (possibly non-Markovian) decays. The final column 'prob' specifies the probability of a random fork and the probabilities of a random fork from a compartment must sum to 1.0. Blank cells in the 'mean' and 'dispersion' column indicate that the cell is copied from the cell above it (i.e., the random fork terminates a decay). Blank cells in the 'prob' column default to 1.0 entries, i.e., to non-random decays. 
<br> The program searches for all filenames (4 filenames, here) starting with the substring **1_add_decay_to_**. Each of the file provides the decays for a different simulation.
3. **2_add_infect_to_** specifies the possible infections for the compartments.
<br> The input is a DataFrame with 4 columns: 'from', 'to', 'infecting', 'R_0'. As one might expect, they specify an infection as a transit by where it starts and ends, and the infecting individual's compartment. The 'R_0' gives the average count of secondary infections from a single infecting individual in a completely susceptible population.
<br> The program searches for all filenames (2 filenames, here) starting with the substring **2_add_infect_to_**. Each of the file provides the infections for a different simulation.
4. **3_add_watchers_1.csv** terminates the algorithm on the threshold day with 30 or more new cases.
**run_ui_threshold.py** has an optional input '-w'. Without a watcher, the algorithm runs to epidemic completion, but with S(0)=1.0e+06 initial susceptibles, e.g., the non-Markovian computations slow the exact Monte Carlo simulation. Although the branching process approximation remains feasible, each of the 24 exact simulations took about 1 hour, because the watcher terminated the algorithm on the day with 30 or more new cases. To explain the watcher input **3_add_watchers_1.csv**, a general explanation of watchers is necessary first.

**Epidemic Watchers**

After the heap-dictionary pops, a user can have any number of watchers examine the present compartment counts and record them. Watchers have 3 types: (1) passive; (2) active; and (3) inactive. The algorithm automatically creates 2 passive watchers, which record compartment counts at certain times: (1) 'all_times', which records at epidemic time 0 and all times after 0 when the heap-dictionary pops; and (2) 'integer_times', which records at epidemic time 0 and and each integer time after. The user can create active watchers to monitor the passive watchers for certain types of conditions, e.g., when a day has 30 or more new cases. An active watcher becomes inactive after its condition is met. The algorithm terminates when every active watcher has been inactivated. It continues to epidemic completion if all watchers are passive.

4. **3_add_watchers_1.csv** terminates the algorithm on the threshold day with 30 or more new cases.
<br>The input is a DataFrame with 5 columns, 'passive_watcher', 'is_difference', 'compartments', 'eta_threshold', and 'outfile'. The DataFrame specifies the active watchers for the simulation.
   1. The column 'passive_watcher' specifies either 'all_times' and 'integer_times', the passive watcher that the active watcher uses.
   2. The column 'is_difference' column specifies either 1 or 0, corresponding to a condition using the difference of compartment counts or just the compartment counts themselves.
   3. The column 'compartment' specifies which compartment counts to monitor, and it can also specify a sum of compartment counts with a '+', e.g., 'i+r' specifies the sum of infectious and recovered compartments.
   4. The column 'eta_threshold' specifies the threshold in the condition, e.g., '29' specifies 30 or more counts.
   5. The column 'outfile' gives the filename of the output CSV file.
  
Putting it all together, the example **3_add_watchers_1.csv** specifies that an active watcher terminate the algorithm when the difference ('is_difference'=1) of sum of compartment counts ('compartments=i+r) at integer times ('watcher'=integer_times) exceeds 29 ('eta_threshold'=29). It specifies the output filename 'pmfs.csv' ('outfile'=pmfs).

The Makefile **Executable/Threshold/run_ui_threshold_make.py** specifies the output directory **Output/Threshold/**, so the output of the active watcher is in the file **Output/Threshold/pmfs.csv**. 

### Code Verification ###

1. The simulation of SIR and SEIR epidemics from a single infectious individual matched analytic formulas for extinction probabilities in branching processes.
2. A compartmental model was initialized with a single I that infects by making S transit directly to R. The distribution of Rs at the end of the epidemic matched Z, the individual reproduction number.
