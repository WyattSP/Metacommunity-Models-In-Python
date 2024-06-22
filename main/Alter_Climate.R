# Reduce Points in a Spectra
# Function to smooth out frequencies
remove.timeSeries.variance <- function(timeseries, lfreq, hfreq, plot.spec = T){
  climate_array = Tmax_4kyr[50,43,]
  lfreq = 0.4
  hfreq = 0.5
  # retrieve current t-series length
  Spectra = timeseries
  n = length(Spectra)
  # check that time series is even length
  if(n%%2 != 0){
    Spectra = Spectra[-1]
    n = length(Spectra)
  } else("Spectra is even")
  # Set frequency interval
  lfreq = lfreq
  hfreq = hfreq
  #############################
  #Extract original properties#
  #############################
  # FFT
  s_fft <- as.complex(fft(Spectra))
  # Find Mean; first value in the series
  mean = abs(s_fft[1])
  # Define the Nyquist frequency
  nyq <- sqrt(Im(s_fft[(n/2) + 1])^2 + Re(s_fft[(n/2) + 1])^2)
  # Define Phases of original series
  phase <- atan2(Im(s_fft), Re(s_fft))
  # Discard first point, nyquist, and second half of series
  s_fft_l = s_fft[2:(n/2)]
  # Define Amplitudes of original series; left side only
  amp <- sqrt(Im(s_fft_l)^2 + Re(s_fft_l)^2)
  # Find Power of original series
  powr <- Re(s_fft_l)^2
  # Portion of imaginary
  phsr <- Im(s_fft_l)
  # Determine frequencies in the original series
  freq <- sapply(1:length(s_fft_l), FUN = function(i) {i/n})
  ###################################################
  #  Remove Power within a Defined Frequency Window #
  ###################################################
  # Set new frequency with additional points
  # Get index values for frequency reset
  zero_freq_index <- which(freq >= lfreq & freq <= hfreq)
  # Change frequencies to zero; phases don't matter
  amp[as.vector(zero_freq_index)] <- 0
  phase[as.vector(zero_freq_index)] <- 0
  #######################
  #      Recompile      #
  #######################
  # Recompile full time series in the frequency domain
  # Need mean, left side of spectrum, nyquist frequency, and right side of spectrum
  # NOTE: Nyquist frequency must be the same as the initial time series to not add in new cycles below 1/f.
  newAseries <- c(mean, amp, nyq, rev(amp))
  newPseries <- c(0, phase[2:((n/2)-1)], 0, phase[((n/2)+1):(length(phase))])
  # NOTE: May want to scale the mean by the ratio of the summed old and new amplitudes
  # Convert to complex
  newSeries = sapply(1:length(newAseries), FUN = function(j) complex(real = newAseries[j]*cos(newPseries[j]), imaginary = (newAseries[j] * sin(newPseries[j]))))
  #######################
  #    Check Lengths    #
  #######################
  L <- length(newSeries)
  if(L != n){
    print("Old and new series length do not match")
  }
  #######################
  #          FFT        #
  #######################
  X = Re(fft(newSeries, inverse = T))/n
  # Plot if desired
  if(plot.spec == TRUE){
    par(mfrow=c(2,2))
    plot(1:length(timeseries),timeseries,col="black",type="l", xlab = "Observations", ylab = "Unit", main = "Original Time Series")
    plot(1:length(X), X,col="red",type="l", xlab = "Observations", ylab = "Unit", main = "Modified Time Series")
    plot(freq,log(powr),type="l")
    plot(freq,log(amp^2),type="l")
    #spectrum(Spectra)
    #spectrum(new_t)
    par(mfrow=c(1,1))
  }
  return(X)
}

# Time series we are working with
Spectra = co2 # define the spectra we want to extend

test_s <- Tmax_4kyr[50,43,]
lfreq_2 = 0.5 - ((1/length(test_s)) * 2)
out_2 <- remove.timeSeries.variance(test_s, lfreq_2, 0.5, T)

lfreq_10 = 0.5 - ((1/length(test_s)) * 10)
out_10 <- remove.timeSeries.variance(test_s, lfreq_10, 0.5, T)

lfreq_25 = 0.5 - ((1/length(test_s)) * 25)
out_25 <- remove.timeSeries.variance(test_s, lfreq_25, 0.5, T)

lfreq_125 = 0.5 - ((1/length(test_s)) * 125)
out_125 <- remove.timeSeries.variance(test_s, lfreq_125, 0.5, T)

lfreq_250 = 0.5 - ((1/length(test_s)) * 250)
out_250 <- remove.timeSeries.variance(test_s, lfreq_250, 0.5, T)

par(mfrow = c(1,5))
spectrum(out_2)
spectrum(out_10)
spectrum(out_25)
spectrum(out_125)
spectrum(out_250)
par(mfrow = c(1,1))
