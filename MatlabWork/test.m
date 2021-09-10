normal_ := plot::Arrow3d(
[0, 0, 0],
[sin(2*a), sin(a)*cos(2*a), cos(a)*cos(2*a)],
a = 0..2*PI):
circle := plot::Circle3d(1, [0, 0, 0], normal_::To,
a = 0..2*PI, Filled):
plot(normal_, circle)
