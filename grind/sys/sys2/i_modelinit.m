function i_modelinit
ax('-defaults','1');
out('-defaults');
disp(' ');
disp('The model is ready for analysis, you can enter GRIND commands now');
disp('Type <a href="matlab: commands">commands</a> to get an overview of all GRIND commands.');
disp('Frequently used commands:');
disp('<a href="matlab: time">time</a>    - run and create time plot');
disp('<a href="matlab: ru">ru</a>      - run and create trajectories in phase plane (or current figure)');
disp('<a href="matlab: null">null</a>    - create null-isoclines in the phase plane');
disp('<a href="matlab: e2p">e2p</a>     - erase figures and re-create time plot and phase plane');
disp('<a href="matlab: paranal">paranal</a> - change one parameter step-by-step and show attractor by simulation');
disp(' ');
