# 2025-slewrate

## Analysis examples

### Time difference between consecutive hits
A simple and satisfactory analysis is the one where we measure the time difference between consecutive hits. For this example, we use [data](https://drive.google.com/file/d/1CjSunUsNMW-gyZt1itdscD9BWwn4I0Vt/view?usp=sharing) collected sending a test-pulse to ALCOR channel 0. An example analysis macro can be found [here](macros/deltat.C).

```
root [0] .L deltat.C
root [1] deltat("alcdaq.fifo_4.root", 0)
```

TDC calibration can be performed with this data. We will see in another example how.
