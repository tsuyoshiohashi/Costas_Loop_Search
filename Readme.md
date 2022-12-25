# Costas Loop Search

# How to use

First,Prepare a file that records the IQ signal of the signal you want to receive.
The reception time should be about 1 second. (If it is long, it will take time to display)
The figure below is an example of creating an IQ file with GNURadio.

![IQ_record](/image/IQ_record.png)

Execute in terminal. Specify the IQ file name as an argument.
~~~
$ python costas_loop_search.py iqfile_name
~~~
After displaying the spectrum for signal confirmation, the symbol-synchronized constellation is displayed.

![PSD](/image/psd.png)
![SymbolSync](/image/after_symbolsync.png)

Although it is circular, it is because the frequency is slightly shifted.
After that, carrier recovery is performed by the Costas loop. 
There are two parameters, alpha and beta, but each of them is changed in 5 ways, calculated, and the results are displayed. It will take some time.
The first plot is the frequency offset. Check the speed until synchronization and the magnitude of jitter at steady state.
The first 3000 symbols are shown for clarity.

![freq_log](/image/freq_log.png)

The next plot similarly plots the constellation with varying parameters.
Since it is BPSK here, two groups of points appear.
The first symbols that are not locked in are removed for clarity.
See if the two clusters of points are clearly separated.

![constellation](/image/constellation.png)

Then select and enter the best parameters from the two plots.
Also, enter the number of symbols until synchronization as an integer.
It is an eye measurement. If anyone knows a better way, please let me know.
The frequency offset and constellation demodulated by the input parameters are displayed.

![freq_log2](/image/freq_log2.png)

![constellation2](/image/constellation2.png)

There is a list of alpha and beta candidates in the program if you want to search by more parameters.
If you rewrite this, it will plot with the new parameters.
When determining the parameters for an actual receiver, it is necessary to check several IQ files.
This is because the frequency offset will vary and fluctuate.

# background

I'm making a BPSK receiver with SDR and Python, but I didn't know how to find good parameters for the Costas loop.
Therefore, I decided to plot the lock-in time and constellation for various parameters and select a value that seems good.

The PySDR documentation (https://pysdr.org/content/rds.html) was very helpful. appreciate.