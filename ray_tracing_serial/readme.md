## Build
Use make to compile

## Run
run like:
	./ray_serial 100000 1000

the first argument is the number of rays, the second is window width

to plot the result, run (using python 3)

python plot.py 

Optionally, you can set the varable create_animation in ray_tracing.c to 1, to create an animation where the light source is moving. Variable number_of_frames and L_v can be modified to affect the animation result. To view the animation, run (using python 3)
python plot_animation.py

