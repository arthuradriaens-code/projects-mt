def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degrees(SomethingInRads):
    return SomethingInRads*180/np.pi
def delta_taccent(theta,deltaz,n):
    v = c/n
    return ((np.cos(theta)*deltaz)/v)*(10**9)


Epsindexofrefractionrange = np.linspace(1.6,1.9,5000)
EpsRelativeAccuracy = np.zeros(len(xcoordinates))
EpsBalloonAngle = np.zeros(len(xcoordinates))

Epstraveltimes = []
Epspaths = []
Epstimes = []
Epsdistances = []

ice = medium.greenland_simple()
prop = radioproparaytracing.radiopropa_ray_tracing(ice, attenuation_model='GL1',config=configh)
for detector in Detectors:
    start_point = Balloon
    final_point = detector
    prop.set_start_and_end_point(start_point, final_point)
    prop.find_solutions()
    SolNumber = prop.get_number_of_solutions()
    for Sol in range(SolNumber):
        Epspaths.append(prop.get_path(Sol))
        Epstimes.append(prop.get_travel_time(Sol))
        x = np.linspace(Epspaths[-1][0,0],start_point[0],1000)
        xlen = x[-1]-x[0]
        z = np.linspace(Epspaths[-1][0,2],start_point[2],1000)
        zlen = z[-1]-z[0]
        diagonallen = np.sqrt(xlen*xlen + zlen*zlen) #(m)
        Epstraveltime = Epstimes[-1]/units.ns + (diagonallen/c)*(10**9) #ns
        Epstraveltimes.append(Epstraveltime)


Epsdifferences = np.zeros(len(Epsindexofrefractionrange))

for number,n in enumerate(Epsindexofrefractionrange):
    #find plane wave
    Epsthetas = np.linspace(0,0.9,1000)
    NumberOfDetectors = len(Detectors)
    Epsdelta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))
    Epsdelta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    Epscorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    Epsnormedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    Epssummedcorrelation = np.zeros(1000)

    for i in range(NumberOfDetectors):
        for j in range(NumberOfDetectors):
            if i < j:
                Epsdelta_t[i][j] = np.abs(traveltimes[i] - traveltimes[j])
                Epsdeltaz = np.linalg.norm(Detectors[i]-Detectors[j])
                Epsposition = (Detectors[0] + Detectors[1])/2
                Epsdelta_taccenten[i][j] = delta_taccent(Epsthetas,np.abs(Epsdeltaz),n)
                Epscorrelation[i][j] = np.abs(Epsdelta_t[i][j] - Epsdelta_taccenten[i][j])
                Epsnormedcorrelation[i][j] = Epscorrelation[i][j]/np.trapz(Epscorrelation[i][j],Epsthetas)

                Epssummedcorrelation += Epsnormedcorrelation[i][j]

    Epsangle_index = np.where(Epssummedcorrelation == Epssummedcorrelation.min())
    Epsangle = Epsthetas[Epsangle_index] #zenith ofc
    b_ballon = MiddleOfDetectors[2]
    a_ballon = (Balloon[2]-b_ballon)/Balloon[0]
    a_planewave = np.tan(np.pi/2-angle)
    b_planewave = MiddleOfDetectors[2]

    angle_snell = np.arcsin(np.sin(angle)*1.27) 
    a_snell = np.tan(np.pi/2-angle_snell)
    b_snell = -1*a_snell*(-1*b_planewave/a_planewave)
    XopBallonHoogte = (Balloon[2] - b_snell)/a_snell
    verschil = XopBallonHoogte - Balloon[0]

    Epsdifferences[number] = np.abs(verschil)

Epsn_index = np.where(Epsdifferences == Epsdifferences.min())
Epsn_fit = Epsindexofrefractionrange[Epsn_index]
if len(n_fit) > 1:
    Epsn_fit = Epsn_fit[0]
Epsdirect_angle = np.pi/2 - np.arctan(a_ballon)
n_actual = ice.get_index_of_refraction(MiddleOfDetectors)

EpsBalloonAngle[s] = degrees(Epsdirect_angle)
print("Epsilon: {}".format(100*(Epsn_fit - n_actual)/n_actual))
