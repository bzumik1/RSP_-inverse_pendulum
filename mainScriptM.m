sys = tf(1 / (s+1))

bode(sys)
nyquist(sys)

pzplot(sys)
step(sys)
impulse(sys)


