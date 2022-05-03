### Results:

- Regarding Kalman Filter estimators, if some states are not measured directly by the sensors, the filter will not be able to estimate the hidden states accurately.
- In this project, I implemented the Unscented Kalman Filter estimator to estimate Radar and Laser measurements. In the case of using only a Laser sensor, the target object velocities cannot be measured directly. So, it causes to increase in the uncertainty of velocity measurements.
- In the case of using only a Radar sensor,  the target object positions cannot be measured directly. So, it causes to increase in the uncertainty of position measurements. 
- As a result of these, fusing Laser and Radar sensors is the best way to estimate states.