
# Description & Usage
- The GitHub respiratory contains the main (Game_of_Life.cpp) and Visual Studio project file. The Game_of_Life.cpp contains initialization parameters in the very beginning (e.g. a Boolean variable for periodicity, an integer variable for the grid size (note: this is the grid size for each process and not the total lattice size – this can be adjusted to roughly get the desired overall grid size), and an integer variable for the total number of iterations. When executed, this will create files with iteration numbers and ids for each process in the debug folder. Merge matrices.py is the main post-processing script which merges files for each iteration into a new file called “iteration_%(iteration_number).txt” as well as outputting the image files “iteration_%(iteration_number).jpg” (note: the user should input the number of iterations, number of rows and columns of the processes – an example is also provided in the GitHub respiratory. Also, all the iteration files must be in the same folder). An animation script animation.py is also provided which uses the jpg image outputs from Mergematrices.py to create a Gif video (example Gif outputs for periodic and non-periodic boundary conditions are provided at GitHub).

![Alt Text](Parrellel-Game-of-Life/Example animation of a 30 by 30 non-periodic lattice (4 processes 15 by 15 each -- 800 iterations)/animation.gif)
