import matplotlib.pyplot as plt

with open('45+') as f:
    lines = f.readlines()
    z_45P = [float(line.split()[0]) for line in lines]
    s_45P = [float(line.split()[9]) for line in lines]
with open('45X') as f:
    lines = f.readlines()
    z_45X = [float(line.split()[0]) for line in lines]
    s_45X = [float(line.split()[9]) for line in lines]
with open('10+') as f:
    lines = f.readlines()
    z_10P = [float(line.split()[0]) for line in lines]
    s_10P = [float(line.split()[9]) for line in lines]
with open('10X') as f:
    lines = f.readlines()
    z_10X = [float(line.split()[0]) for line in lines]
    s_10X = [float(line.split()[9]) for line in lines]
with open('15+') as f:
    lines = f.readlines()
    z_15P = [float(line.split()[0]) for line in lines]
    s_15P = [float(line.split()[9]) for line in lines]
with open('15X') as f:
    lines = f.readlines()
    z_15X = [float(line.split()[0]) for line in lines]
    s_15X = [float(line.split()[9]) for line in lines]
with open('20+') as f:
    lines = f.readlines()
    z_20P = [float(line.split()[0]) for line in lines]
    s_20P = [float(line.split()[9]) for line in lines]
with open('20X') as f:
    lines = f.readlines()
    z_20X = [float(line.split()[0]) for line in lines]
    s_20X = [float(line.split()[9]) for line in lines]
with open('25+') as f:
    lines = f.readlines()
    z_25P = [float(line.split()[0]) for line in lines]
    s_25P = [float(line.split()[9]) for line in lines]
with open('25X') as f:
    lines = f.readlines()
    z_25X = [float(line.split()[0]) for line in lines]
    s_25X = [float(line.split()[9]) for line in lines]
with open('30+') as f:
    lines = f.readlines()
    z_30P = [float(line.split()[0]) for line in lines]
    s_30P = [float(line.split()[9]) for line in lines]
with open('30X') as f:
    lines = f.readlines()
    z_30X = [float(line.split()[0]) for line in lines]
    s_30X = [float(line.split()[9]) for line in lines]
with open('35+') as f:
    lines = f.readlines()
    z_35P = [float(line.split()[0]) for line in lines]
    s_35P = [float(line.split()[9]) for line in lines]
with open('35X') as f:
    lines = f.readlines()
    z_35X = [float(line.split()[0]) for line in lines]
    s_35X = [float(line.split()[9]) for line in lines]
with open('40+') as f:
    lines = f.readlines()
    z_40P = [float(line.split()[0]) for line in lines]
    s_40P = [float(line.split()[9]) for line in lines]
with open('40X') as f:
    lines = f.readlines()
    z_40X = [float(line.split()[0]) for line in lines]
    s_40X = [float(line.split()[9]) for line in lines]
with open('50+') as f:
    lines = f.readlines()
    z_50P = [float(line.split()[0]) for line in lines]
    s_50P = [float(line.split()[9]) for line in lines]
with open('50X') as f:
    lines = f.readlines()
    z_50X = [float(line.split()[0]) for line in lines]
    s_50X = [float(line.split()[9]) for line in lines]
with open('60+') as f:
    lines = f.readlines()
    z_60P = [float(line.split()[0]) for line in lines]
    s_60P = [float(line.split()[9]) for line in lines]
with open('60X') as f:
    lines = f.readlines()
    z_60X = [float(line.split()[0]) for line in lines]
    s_60X = [float(line.split()[9]) for line in lines]
with open('70+') as f:
    lines = f.readlines()
    z_70P = [float(line.split()[0]) for line in lines]
    s_70P = [float(line.split()[9]) for line in lines]
with open('70X') as f:
    lines = f.readlines()
    z_70X = [float(line.split()[0]) for line in lines]
    s_70X = [float(line.split()[9]) for line in lines]
with open('80+') as f:
    lines = f.readlines()
    z_80P = [float(line.split()[0]) for line in lines]
    s_80P = [float(line.split()[9]) for line in lines]
with open('80X') as f:
    lines = f.readlines()
    z_80X = [float(line.split()[0]) for line in lines]
    s_80X = [float(line.split()[9]) for line in lines]

plt.figure(layout="constrained")


ax1 = plt.subplot(349)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_50P,s_50P,color='green')
plt.scatter(z_50X,s_50X,color='red')
plt.legend(['Sphere size= 50m'])
plt.ylabel("relative computational time")
plt.xlabel("stepsize (°)")

ax2 = plt.subplot(3,4,10)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_60P,s_60P,color='green')
plt.scatter(z_60X,s_60X,color='red')
plt.legend(['Sphere size= 60m'])
plt.tick_params('y', labelleft=False)
plt.xlabel("stepsize (°)")

ax3 = plt.subplot(3,4,11)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_70P,s_70P,color='green')
plt.scatter(z_70X,s_70X,color='red')
plt.legend(['Sphere size= 70m'])
plt.tick_params('y', labelleft=False)
plt.xlabel("stepsize (°)")

ax4 = plt.subplot(3,4,12)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_80P,s_80P,color='green')
plt.scatter(z_80X,s_80X,color='red')
plt.legend(['Sphere size= 80m'])
plt.tick_params('y', labelleft=False)
plt.xlabel("stepsize (°)")

plt.subplot(341,sharex=ax1)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_10P,s_10P,color='green')
plt.scatter(z_10X,s_10X,color='red')
plt.legend(['Sphere size= 10m'])
plt.tick_params('x', labelbottom=False)
plt.ylabel("relative computational time")

plt.subplot(342)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_15P,s_15P,color='green')
plt.scatter(z_15X,s_15X,color='red')
plt.legend(['Sphere size= 15m'])
plt.tick_params('x', labelbottom=False)
plt.tick_params('y', labelleft=False)

plt.subplot(343)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_20P,s_20P,color='green')
plt.scatter(z_20X,s_20X,color='red')
plt.legend(['Sphere size= 20m'])
plt.tick_params('x', labelbottom=False)
plt.tick_params('y', labelleft=False)

plt.subplot(344)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_25P,s_25P,color='green')
plt.scatter(z_25X,s_25X,color='red')
plt.legend(['Sphere size= 25m'])
plt.tick_params('x', labelbottom=False)
plt.tick_params('y', labelleft=False)

plt.subplot(345,sharex=ax1)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_30P,s_30P,color='green')
plt.scatter(z_30X,s_30X,color='red')
plt.legend(['Sphere size= 30m'])
plt.tick_params('x', labelbottom=False)
plt.ylabel("relative computational time")

plt.subplot(346)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_35P,s_35P,color='green')
plt.scatter(z_35X,s_35X,color='red')
plt.legend(['Sphere size= 35m'])
plt.tick_params('x', labelbottom=False)
plt.tick_params('y', labelleft=False)

plt.subplot(348)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_45P,s_45P,color='green')
plt.scatter(z_45X,s_45X,color='red')
plt.legend(['Sphere size= 45m'])
plt.tick_params('x', labelbottom=False)
plt.tick_params('y', labelleft=False)


ax = plt.subplot(347)
plt.xlim(0,1.2)
plt.ylim(0,2.4)
plt.scatter(z_40P,s_40P,color='green')
plt.scatter(z_40X,s_40X,color='red')
plt.legend(['Sphere size= 40m'])
plt.tick_params('x', labelbottom=False)
plt.tick_params('y', labelleft=False)


plt.show()
