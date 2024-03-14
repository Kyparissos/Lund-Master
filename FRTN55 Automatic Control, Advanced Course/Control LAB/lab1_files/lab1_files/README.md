This forlder contains the matlab files for the mass-spring laboratory in the
 Advanced Control course (previoously Multivariable Control).


 * define_simulation_process.m : used in the servo_simulated.mdl by the block defining the model parameters
 * init_servo.m : used by the simulink models to initialize controller and process
 * loop_shaping.m : script for control design. Also used to define the feedback and feed-forward controllers.
 * servo_model.m : define the simulation process of the mass-spring system. Needs to bve run before the loop_shaping.m script.
 * servo_real.mdl : simulink model to interface and control the physical process.
 * servo_simulated.mdl : simulink model to simulate the whole control system.
 * specs.m : used in the servo_simulated.mdl by the block "Check specifications". It plots the results of the simulaion against the time-domain specifications.
 * specs_real.m : used in the servo_real.mdl by the block "Check specifications". It plots the results of the test against the time-domain specifications.

