using PyPlot

a=2; b=4

x = range(-10, stop=10, length=5000)
y = range(-10, stop=10, length=5000)

z=@. (x^4)/a^4+(y/b)'^4-(x'^2*y^2)/(a*b)
surf(x, y, z,alpha = 0.9, cmap="inferno")
PyPlot.display_figs()
     