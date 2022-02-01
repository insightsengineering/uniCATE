## uniCATE 0.3.0 (2022-02-01)

* `uniCATE` is going public!

## uniCATE 0.2.3 (2022-01-31)

* Adding GitHub Actions for CI and CD.
* Made largely stylistic updates to internal functions to pass linting and style
  checks.
* Set the license to Apache 2.0


## uniCATE 0.2.2 (2022-01-31)

* `failure` argument in `sunicate()` has been changed to `event`.
* Change to `sunicate()` behaviour: The median of the `data` argument's
  `relative_time` variable is now set as the default `time_cutoff` when
  `time_cutoff` is `NULL`.
* Change to internal `sunicate()` behaviour: When transforming wide data to a
  long format, no more than five unique relative times are reported for each
  observation. These relative times correspond to the quintiles of the relative
  times between the earliest relative time, and the minimum of each observations
  relative time and the time cutoff.

## uniCATE 0.2.1 (2022-01-21)

* Documentation of internal and exported functions has been cleaned.
* New package introduction in README.
* Co-authors added to DESCRIPTION
* Package description updated in DESCRIPTION

## uniCATE 0.2.0 (2022-01-18)

* Added `sunicate()`, the `unicate()` counterpart for right-censored
  time-to-event outcomes.
* Minor touch-ups to miscellaneous docs. 
