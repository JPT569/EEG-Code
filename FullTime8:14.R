
cat("Power + Autocorr + Entropy \n")

# Load required libraries
library(edfReader)
library(ggplot2)
library(data.table)
library(signal)
library(multitaper)
library(e1071)
library(tidyr)

setwd("/Users/jacobtiller/Desktop/EEG Biomarker Work")
theFile <- "SLC6A1 WCM07 10.6h 9896.edf"

edf_header <- readEdfHeader(theFile)
edf_signals <- readEdfSignals(edf_header)
sampleRate <- edf_signals$C3$sRate

recording_start_exact <- as.POSIXct(edf_header$startTime, format="%Y-%m-%d %H:%M:%S")
cat(sprintf("Recording started: %s\n", format(recording_start_exact, "%H:%M:%S")))

# Channel mapping
find_channel_name <- function(target_name, available_names) {
  if(target_name %in% available_names) return(target_name)
  
  variants <- c(target_name, toupper(target_name), tolower(target_name),
                paste0(toupper(substr(target_name, 1, 1)), tolower(substr(target_name, 2, nchar(target_name)))))
  
  if(target_name == "T3") variants <- c(variants, "T7")
  if(target_name == "T4") variants <- c(variants, "T8")
  if(target_name == "T5") variants <- c(variants, "P7")
  if(target_name == "T6") variants <- c(variants, "P8")
  if(target_name == "T7") variants <- c(variants, "T3")
  if(target_name == "T8") variants <- c(variants, "T4")
  
  for(variant in variants) {
    if(variant %in% available_names) return(variant)
  }
  return(NULL)
}

# Define channel pairs
channel_pairs_raw <- list(
  c("FP1", "F3"), c("F3", "C3"), c("C3", "P3"), c("P3", "O1"),
  c("FP2", "F4"), c("F4", "C4"), c("C4", "P4"), c("P4", "O2"),
  c("FP1", "F7"), c("F7", "T3"), c("T3", "T5"), c("T5", "O1"),
  c("FP2", "F8"), c("F8", "T4"), c("T4", "T6"), c("T6", "O2"),
  c("FZ", "CZ"), c("CZ", "PZ")
)

channel_pairs <- list()
for(i in 1:length(channel_pairs_raw)) {
  pair_raw <- channel_pairs_raw[[i]]
  ch1_mapped <- find_channel_name(pair_raw[1], names(edf_signals))
  ch2_mapped <- find_channel_name(pair_raw[2], names(edf_signals))
  
  if(!is.null(ch1_mapped) && !is.null(ch2_mapped)) {
    channel_pairs[[length(channel_pairs) + 1]] <- c(ch1_mapped, ch2_mapped)
  }
}

channel_pair_names <- sapply(channel_pairs, function(x) paste(x[1], x[2], sep="-"))

# Original + New seizure zones
original_seizure_times <- data.frame(
  label = c("18:23:57-05", "18:24:14-17", "18:24:20", "18:24:26-27", "18:24:32-38", "18:24:40-47", 
            "18:24:53", "18:25:49-50", "18:26:19-26", "18:26:34-35", "18:26:41-47"),
  start_time_str = c("18:23:57", "18:24:14", "18:24:20", "18:24:26", "18:24:32", "18:24:40", 
                     "18:24:53", "18:25:49", "18:26:19", "18:26:34", "18:26:41"),
  end_time_str = c("18:24:06", "18:24:18", "18:24:21", "18:24:28", "18:24:39", "18:24:48", 
                   "18:24:54", "18:25:51", "18:26:27", "18:26:36", "18:26:48"),
  color_type = c("green", "green", "green", "green", "green", "green", 
                 "yellow", "green", "green", "green", "green"),
  stringsAsFactors = FALSE
)

new_seizure_times <- data.frame(
  label = c("18:26:56-18:27:01", "18:27:12-18:27:19", "18:27:27-18:27:28", "18:27:30", 
            "18:27:42-18:27:45", "18:27:54-18:28:03", "18:28:04-18:28:05", "18:28:12-18:28:13",
            "18:28:22", "18:28:46-18:28:47", "18:29:15-18:29:20", "18:29:23-18:29:24",
            "18:29:34", "18:29:43-18:29:48", "18:29:53-18:29:55", "18:30:11-18:30:15"),
  start_time_str = c("18:26:56", "18:27:12", "18:27:27", "18:27:30", 
                     "18:27:42", "18:27:54", "18:28:04", "18:28:12",
                     "18:28:22", "18:28:46", "18:29:15", "18:29:23",
                     "18:29:34", "18:29:43", "18:29:53", "18:30:11"),
  end_time_str = c("18:27:01", "18:27:19", "18:27:28", "18:27:30", 
                   "18:27:45", "18:28:03", "18:28:05", "18:28:13",
                   "18:28:22", "18:28:47", "18:29:20", "18:29:24",
                   "18:29:34", "18:29:48", "18:29:55", "18:30:15"),
  color_type = c("green", "green", "yellow", "yellow", 
                 "green", "green", "yellow", "yellow",
                 "yellow", "yellow", "yellow", "yellow",
                 "green", "green", "green", "green"),
  stringsAsFactors = FALSE
)

# Combine and exclude specified periods
combined_seizure_times <- rbind(original_seizure_times, new_seizure_times)
excluded_labels <- c("18:28:46-18:28:47", "18:29:15-18:29:20")

# Time conversion
convert_to_offset <- function(time_str) {
  target_time <- as.POSIXct(paste(as.Date(recording_start_exact), time_str), format="%Y-%m-%d %H:%M:%S")
  offset <- as.numeric(difftime(target_time, recording_start_exact, units = "secs"))
  return(offset)
}

# Create corrected zones (excluding specified periods)
all_zones <- data.frame()
for(i in 1:nrow(combined_seizure_times)) {
  start_offset <- convert_to_offset(combined_seizure_times$start_time_str[i])
  end_offset <- convert_to_offset(combined_seizure_times$end_time_str[i])
  
  all_zones <- rbind(all_zones, data.frame(
    label = combined_seizure_times$label[i],
    start_time = start_offset,
    end_time = end_offset,
    start_str = combined_seizure_times$start_time_str[i],
    end_str = combined_seizure_times$end_time_str[i],
    color_type = combined_seizure_times$color_type[i]
  ))
}

# Remove excluded periods entirely
correct_zones <- all_zones[!all_zones$label %in% excluded_labels, ]

cat(sprintf("Corrected zones: %d (excluded %s)\n", nrow(correct_zones), paste(excluded_labels, collapse=", ")))

analyze_eeg_final <- function(theSignal, sampling_rate, power_low, power_high, powerWindowLengthSeconds, correlationWindowLengthSeconds, stepSizeSeconds) {
  n_samples <- length(theSignal)
  power_window_length <- sampling_rate * powerWindowLengthSeconds
  step_size <- round(stepSizeSeconds * sampling_rate)
  n_windows <- floor((n_samples - power_window_length) / step_size) + 1
  
  power_values <- numeric(n_windows)
  full_entropy <- numeric(n_windows)
  
  for (i in 1:n_windows) {
    start_idx <- (i - 1) * step_size + 1
    end_idx <- start_idx + power_window_length - 1
    segment <- theSignal[start_idx:end_idx]
    
    # Power spectral density
    psd <- spectrum(segment, plot = FALSE, method = "pgram", spans = c(5))
    freqs <- psd$freq * sampling_rate
    power_values[i] <- sum(psd$spec[freqs >= power_low & freqs <= power_high])
    
    # Full spectrum entropy (1-30 Hz)
    freq_mask_full <- freqs >= 1 & freqs <= 30
    if (sum(freq_mask_full) > 1) {
      full_psd <- psd$spec[freq_mask_full]
      full_psd_norm <- full_psd / sum(full_psd)
      full_psd_norm[full_psd_norm <= 0] <- 1e-10
      full_entropy[i] <- -sum(full_psd_norm * log2(full_psd_norm))
    } else {
      full_entropy[i] <- NA
    }
  }
  
  # Autocorrelation
  autocorr_window_length <- correlationWindowLengthSeconds * sampling_rate
  lag_min <- 1 / power_high
  lag_max <- 1 / power_low
  lags <- round(seq(lag_min, lag_max, length.out = 10) * sampling_rate)
  
  max_autocorrelation <- numeric(n_windows)
  for (i in 1:n_windows) {
    start_idx <- (i - 1) * step_size + 1
    end_idx <- start_idx + autocorr_window_length - 1
    if (end_idx > n_samples) end_idx <- n_samples
    
    segment <- theSignal[start_idx:end_idx]
    autocorr <- sapply(lags, function(lag) {
      if(lag < length(segment)) {
        cor(segment[1:(length(segment) - lag)], segment[(lag + 1):length(segment)], use="complete.obs")
      } else {
        0
      }
    })
    max_autocorrelation[i] <- max(autocorr, na.rm = TRUE)
  }
  
  time_vector_windows <- seq(0, (n_windows - 1) * stepSizeSeconds, by = stepSizeSeconds)
  
  # Apply shifts
  shiftSignal <- function(x, theShift) {
    return(c(rep(0, theShift), x[1:(length(x) - theShift)]))
  }
  
  power_values <- shiftSignal(power_values, floor((powerWindowLengthSeconds / stepSizeSeconds)/2))
  max_autocorrelation <- shiftSignal(max_autocorrelation, floor(correlationWindowLengthSeconds / stepSizeSeconds / 2))
  full_entropy <- shiftSignal(full_entropy, floor((powerWindowLengthSeconds / stepSizeSeconds)/2))
  
  result_data <- data.table(
    Time = time_vector_windows,
    Power = power_values,
    Max_Autocorrelation = max_autocorrelation,
    Full_Entropy = full_entropy
  )
  
  return(result_data)
}

cat("\nThis is all channels\n")

nyquist_freq <- sampleRate / 2
low <- 1 / nyquist_freq
high <- 70 / nyquist_freq
butter_bandpass <- butter(4, c(low, high), type = "pass")

clinical_channel_data <- list()

for(i in 1:length(channel_pairs)) {
  pair <- channel_pairs[[i]]
  pair_name <- channel_pair_names[i]
  
  signal1_name <- pair[1]
  signal2_name <- pair[2]
  
  if(signal1_name %in% names(edf_signals) && signal2_name %in% names(edf_signals)) {
    signal1 <- edf_signals[[signal1_name]]$signal
    signal2 <- edf_signals[[signal2_name]]$signal
    
    diff_signal <- signal1 - signal2
    diff_signal[is.na(diff_signal)] <- 0
    
    diff_signal_filtered <- filtfilt(butter_bandpass, diff_signal)
    diff_signal_filtered_120_mins <- na.omit(diff_signal_filtered[1:(sampleRate*60*120)])
    
    result <- analyze_eeg_final(diff_signal_filtered_120_mins, sampleRate, 2.5, 4, 1, 1, 0.25)
    
    # Apply OPTIMIZED thresholds (from stricterrangesfull)
    opt_power_thresh <- 105932      
    opt_autocorr_thresh <- 0.407    
    opt_entropy_thresh <- 3.16      
    
    result[, opt_baseline_power := (Power > opt_power_thresh)]
    result[, opt_baseline_autocorr := (Max_Autocorrelation > opt_autocorr_thresh)]
    result[, opt_baseline_entropy := (Full_Entropy <= opt_entropy_thresh)]
    
    # OPTIMIZED METHOD: Power + Autocorr + Entropy (ALL 3 required)
    result[, methodA_opt_current := opt_baseline_power & opt_baseline_autocorr & opt_baseline_entropy]
    
    clinical_channel_data[[pair_name]] <- result
    
    cat(sprintf("Processed %s: %d time points\n", pair_name, nrow(result)))
  }
}

cat(sprintf("Processed %d channel pairs\n", length(clinical_channel_data)))

plot_start <- convert_to_offset("18:23:33")
plot_end <- convert_to_offset("18:30:26")

cat(sprintf("\nAnalysis window: 18:23:33 to 18:30:26\n"))

# Enhanced evaluation function
evaluate_optimized_method <- function() {
  total_hits <- 0
  hits_in_zones <- 0
  
  for(pair_name in names(clinical_channel_data)) {
    pair_data <- clinical_channel_data[[pair_name]]
    analysis_window <- pair_data[Time >= plot_start & Time <= plot_end]
    
    method_hits <- analysis_window[methodA_opt_current == TRUE]
    
    if(nrow(method_hits) > 0) {
      for(i in 1:nrow(method_hits)) {
        hit_time <- method_hits$Time[i]
        total_hits <- total_hits + 1
        
        # Check if in corrected seizure zone (±2s buffer)
        is_in_zone <- FALSE
        for(j in 1:nrow(correct_zones)) {
          zone_start_extended <- correct_zones$start_time[j] - 2
          zone_end_extended <- correct_zones$end_time[j] + 2
          
          if(hit_time >= zone_start_extended & hit_time <= zone_end_extended) {
            is_in_zone <- TRUE
            break
          }
        }
        
        if(is_in_zone) {
          hits_in_zones <- hits_in_zones + 1
        }
      }
    }
  }
  
  # Calculate sensitivity
  seizure_windows_covered <- 0
  total_seizure_windows <- 0
  
  # Create time grid at 0.25s resolution
  time_points <- seq(plot_start, plot_end, by = 0.25)
  
  # First identify all seizure time points
  for(time_point in time_points) {
    is_seizure_time <- FALSE
    for(j in 1:nrow(correct_zones)) {
      if(time_point >= correct_zones$start_time[j] && time_point <= correct_zones$end_time[j]) {
        is_seizure_time <- TRUE
        break
      }
    }
    
    if(is_seizure_time) {
      total_seizure_windows <- total_seizure_windows + 1
      
      # Check if this seizure time point is covered by any detection within ±2 seconds
      is_detected <- FALSE
      for(pair_name in names(clinical_channel_data)) {
        pair_data <- clinical_channel_data[[pair_name]]
        nearby_hits <- pair_data[methodA_opt_current == TRUE & 
                                   Time >= plot_start & Time <= plot_end &
                                   abs(Time - time_point) <= 2]
        if(nrow(nearby_hits) > 0) {
          is_detected <- TRUE
          break
        }
      }
      
      if(is_detected) {
        seizure_windows_covered <- seizure_windows_covered + 1
      }
    }
  }
  
  sensitivity <- ifelse(total_seizure_windows > 0, 100 * seizure_windows_covered / total_seizure_windows, 0)
  
  ppv <- ifelse(total_hits > 0, 100 * hits_in_zones / total_hits, 0)
  
  return(list(
    total_hits = total_hits,
    hits_in_zones = hits_in_zones,
    ppv = ppv,
    sensitivity = sensitivity,
    total_seizure_windows = total_seizure_windows,
    covered_seizure_windows = seizure_windows_covered
  ))
}

# Evaluate optimized method
method_result <- evaluate_optimized_method()

cat("\nPerformance\n")
cat(sprintf("Power + Autocorr + Entropy (Optimized):\n"))
cat(sprintf("  PPV: %.1f%%\n", method_result$ppv))
cat(sprintf("  Sensitivity: %.1f%%\n", method_result$sensitivity))
cat(sprintf("  Total hits: %d\n", method_result$total_hits))
cat(sprintf("  Hits in zones: %d\n", method_result$hits_in_zones))

# Collect all detection hits
all_hits <- data.table()

for(pair_name in names(clinical_channel_data)) {
  pair_data <- clinical_channel_data[[pair_name]]
  
  method_hits <- pair_data[methodA_opt_current == TRUE & Time >= plot_start & Time <= plot_end]
  if(nrow(method_hits) > 0) {
    method_hits[, Channel := pair_name]
    all_hits <- rbind(all_hits, method_hits[, .(Time, Channel)])
  }
}

cat(sprintf("Hits in analysis window: %d\n", nrow(all_hits)))

# Assign Y positions
unique_channels <- unique(all_hits$Channel)
if(length(unique_channels) == 0) {
  unique_channels <- names(clinical_channel_data)
}

if(nrow(all_hits) > 0) {
  all_hits[, Y_Position := match(Channel, unique_channels)]
}

# Time setup
time_breaks <- seq(plot_start, plot_end, by = 30)
time_labels <- sapply(time_breaks, function(offset) {
  real_time <- recording_start_exact + offset
  format(real_time, "%H:%M:%S")
})

# Separate zones by color
green_zones <- correct_zones[correct_zones$color_type == "green", ]
yellow_zones <- correct_zones[correct_zones$color_type == "yellow", ]

cat("\n Plot \n")

optimized_plot <- ggplot() +
  # Corrected seizure zones
  geom_rect(data = green_zones, 
            aes(xmin = start_time, xmax = end_time, ymin = 0.5, ymax = length(unique_channels) + 0.5),
            fill = "green", alpha = 0.3) +
  geom_rect(data = yellow_zones, 
            aes(xmin = start_time, xmax = end_time, ymin = 0.5, ymax = length(unique_channels) + 0.5),
            fill = "yellow", alpha = 0.4) +
  
  # Detection hits (blue triangles)
  {if(nrow(all_hits) > 0) {
    geom_point(data = all_hits, 
               aes(x = Time, y = Y_Position), 
               color = "blue", size = 2.0, alpha = 0.8, shape = 17)
  }} +
  
  scale_y_continuous(name = "Channel Pair", 
                     breaks = 1:length(unique_channels), 
                     labels = unique_channels,
                     limits = c(0.5, length(unique_channels) + 0.5)) +
  scale_x_continuous(name = "Time", 
                     breaks = time_breaks, 
                     labels = time_labels,
                     limits = c(plot_start, plot_end)) +
  
  labs(title = "POWER + AUTOCORR + ENTROPY (Optimized Thresholds)",
       subtitle = sprintf("Blue = Detections (%d hits) | PPV = %.1f%% | Sensitivity = %.1f%% | Zones = %d", 
                          nrow(all_hits), 
                          method_result$ppv,
                          method_result$sensitivity,
                          nrow(correct_zones)),
       caption = "Power>105932, Autocorr>0.407, Entropy<=3.16 | ALL 3 required | Excludes 18:29:15-20 and 18:28:46-47") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.minor = element_blank())

# Add zone labels
for(i in 1:nrow(correct_zones)) {
  zone_center <- (correct_zones$start_time[i] + correct_zones$end_time[i]) / 2
  label_color <- ifelse(correct_zones$color_type[i] == "yellow", "orange", "darkgreen")
  optimized_plot <- optimized_plot + annotate("text", x = zone_center, y = length(unique_channels) + 0.2, 
                                              label = correct_zones$start_str[i], color = label_color, size = 2.5, fontface = "bold")
}

# Display the plot
print(optimized_plot)

cat(sprintf("\n FINAL RESULTS:\n"))
cat(sprintf("• Method: Power + Autocorr + Entropy (ALL 3 required)\n"))
cat(sprintf("• PPV: %.1f%% | Sensitivity: %.1f%% | Hits: %d\n", 
            method_result$ppv, method_result$sensitivity, method_result$total_hits))
cat(sprintf("• Thresholds: Power>105932, Autocorr>0.407, Entropy<=3.16\n"))
cat(sprintf("• Zones: %d corrected seizure periods\n", nrow(correct_zones)))
