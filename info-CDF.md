# CDF file

## AIA Standard [ANDI] Variables

- ordinate_values [float32] -> Intensity [array]. These are the detector signal values for each point.
- actual_delay_time [float32] -> Time since injection until the beginning of detector measurement (often referred to as the "dead time" or programmed start delay).
- actual_sampling_interval [float32] -> Time between data points. (Example: A value of 0.5 means the detector recorded a signal every 0.5 seconds).

## Axis reconstruction

- X-axis -> time[n] = actual_delay_time + (n * actual_sampling_interval), where n is the index of point
- X-axis -> ordinate_values

## Units

- Intensity -> Arbitrary
- Retention time -> Seconds
