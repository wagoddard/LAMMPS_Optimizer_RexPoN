
# done in ipython -pylab mode

x = array([10,15,20,25])
y = array([450,1334,3194,5888])

z = poly1d( polyfit(x,y,3) )         # FITTING HERE, z is evaluating the polynimial

az = poly1d( array([0.2, 0, 0, 0]) ) # my approximate fit 0.2 * x^3

cx = linspace(min(nx), max(nx), 100)

plot(x, y, '.', cx, z(cx), '-', cx, az(cx), '--')


