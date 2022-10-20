
    def hybrid_minimizer(self,n_reflections=0):
        """
        Uses RadioPropa to first find all the numerical ray tracing solutions between sphere x1 
        and sphere x2 for a big sphere. After which it uses the scipy optimize.brentq module
        to find the best path in these angle intervals. 
        Tracer does not work for reflective bottoms or secondary creation at the moment
        """

        try:
            X1 = self._X1 * (radiopropa.meter/units.meter)
            X2 = self._X2 * (radiopropa.meter/units.meter)
        except TypeError: 
            self.__logger.error('NoneType: start or endpoint not initialized')
            raise TypeError('NoneType: start or endpoint not initialized')
      
        v = (self._X2 - self._X1)
        u = copy.deepcopy(v)
        u[2] = 0
        theta_direct, phi_direct = hp.cartesian_to_spherical(*v) # zenith and azimuth for the direct linear ray solution (radians)
        cherenkov_angle = np.arccos(1. / self._medium.get_index_of_refraction(self._X1))
        
        ## regions of theta with posible solutions (radians)
        launch_lower = [0]
        launch_upper = [theta_direct + 2*abs(self.delta_theta_direct(dz=self._sphere_sizes[0]))] # below theta_direct no solutions are possible without upward reflections


        detected_rays = []
        detected_theta = []
        results = []
        res_angle = 0.001*units.degree/units.radian

        def get_ray(theta,phi):
            ray_dir = hp.spherical_to_cartesian(theta,phi)
            if ray_dir.shape==(3,1): ray_dir = ray_dir.T[0] #doesn't always give the right shape
            source = radiopropa.Source()
            source.add(radiopropa.SourcePosition(radiopropa.Vector3d(*x1)))
            source.add(radiopropa.SourceDirection(radiopropa.Vector3d(*ray_dir)))
            ray = source.getCandidate()
            return ray

        def shoot_ray(theta):
            ray = get_ray(theta,phi_direct)
            sim.run(ray, True)
            return ray

        def cot(x):
            return 1/np.tan(x)

        def arccot(x):
            return np.arctan(-x) + np.pi/2

        def delta_z(cot_theta):
            theta = arccot(cot_theta)
            ray = shoot_ray(theta)
            ray_endpoint = self.get_path_candidate(ray)[-1]
            return (ray_endpoint-self._X2)[2]

        def delta_z_squared(cot_theta):
            return delta_z(cot_theta)**2
        

        if n_reflections > 0:
            if self.medium.reflection is None:
                self.__logger.error("a solution for {:d} reflection(s) off the bottom reflective layer is requested,"
                                        +"but ice model does not specify a reflective layer".format(n_reflections))
                raise AttributeError("a solution for {:d} reflection(s) off the bottom reflective layer is requested,"
                                        +"but ice model does not specify a reflective layer".format(n_reflections))
            else:
                z_refl = self._medium.reflection
                rho_channel = np.linalg.norm(u)
                if self._X2[2] > self._X1[2]: 
                    z_up = self._X2[2]
                    z_down = self._X1[2]
                else:
                    z_up = self._X1[2]
                    z_down = self._X2[2]
                rho_bottom = (rho_channel * (z_refl - z_down)) / (2*z_refl - z_up - z_down)
                alpha = np.arctan((z_down - z_refl)/rho_bottom)
                ## when reflection on the bottom are allowed, a initial region for theta from 180-alpha to 180 degrees is added
                launch_lower.append(((np.pi/2 + alpha) - 2*abs(self.delta_theta_bottom(dz=self._sphere_sizes[0], z_refl=z_refl) / units.radian)))
                launch_upper.append(np.pi)

        sphere_size = self._sphere_sizes[0] * (radiopropa.meter/units.meter)
        s = 0
        detected_rays = []
        results = []

        ##define module list for simulation
        sim = radiopropa.ModuleList()
        sim.add(radiopropa.PropagationCK(self._ice_model.get_scalar_field(), 1E-8, .001, 1.)) ## add propagation to module list
        for module in self._ice_model.get_modules().values(): 
            sim.add(module)
        sim.add(radiopropa.MaximumTrajectoryLength(self._max_traj_length * (radiopropa.meter/units.meter)))

        ## define observer for detection (channel)            
        obs = radiopropa.Observer()
        obs.setDeactivateOnDetection(True)
        channel = radiopropa.ObserverSurface(radiopropa.Sphere(radiopropa.Vector3d(*X2), sphere_size)) ## when making the radius larger than 2 meters, somethimes three solution times are found
        obs.add(channel)
        sim.add(obs)

        ## define observer for stopping simulation (boundaries)
        obs2 = radiopropa.Observer()
        obs2.setDeactivateOnDetection(True)
        w = (u / np.linalg.norm(u)) * 2*sphere_size
        boundary_behind_channel = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(*(X2 + w)), radiopropa.Vector3d(*w)))
        obs2.add(boundary_behind_channel)
        boundary_above_surface = radiopropa.ObserverSurface(radiopropa.Plane(radiopropa.Vector3d(0, 0, 1*radiopropa.meter), radiopropa.Vector3d(0, 0, 1)))
        obs2.add(boundary_above_surface)
        sim.add(obs2)

        t1 = time.time()
        
        #create total scanning range from the upper and lower thetas of the bundles
        step = self._step_zeniths[s] / units.radian
        theta_scanning_range = np.array([])
        for iL in range(len(launch_lower)):
            new_scanning_range = np.arange(launch_lower[iL], launch_upper[iL]+step, step)
            theta_scanning_range = np.concatenate((theta_scanning_range, new_scanning_range))

        for theta in theta_scanning_range:
            ray_dir = hp.spherical_to_cartesian(theta, phi_direct)
            
            def delta(ray_dir,shower_dir):
                viewing = np.arccos(np.dot(shower_dir, ray_dir)) * units.radian
                return viewing - cherenkov_angle

            if (self._shower_axis is None) or (abs(delta(ray_dir,self._shower_axis)) < self._cut_viewing_angle):
                source = radiopropa.Source()
                source.add(radiopropa.SourcePosition(radiopropa.Vector3d(*X1)))
                source.add(radiopropa.SourceDirection(radiopropa.Vector3d(*ray_dir)))
                sim.setShowProgress(True)
                ray = source.getCandidate()
                sim.run(ray, True)
                
                current_rays = [ray]
                while len(current_rays) > 0:
                    next_rays = []
                    for ray in current_rays:
                        if channel.checkDetection(ray.get()) == radiopropa.DETECTED:
                            detected_rays.append(ray)
                            result = {}
                            if n_reflections == 0:
                                result['reflection']=0
                                result['reflection_case']=1
                            elif self._ice_model.get_modules()["bottom reflection"].get_times_reflectedoff(ray.get()) <= n_reflections: 
                                result['reflection']=self._ice_model.get_modules()["bottom reflection"].get_times_reflectedoff(ray.get())
                                result['reflection_case']=int(np.ceil(theta/np.deg2rad(90)))
                            results.append(result)
                        for secondary in ray.secondaries:
                            next_rays.append(secondary)
                    current_rays = next_rays

        #loop over previous rays to find the upper and lower theta of each bundle of rays
        #uses step, but because step is initialized after this loop this is the previous step size as intented
        if len(detected_rays) > 0:
            launch_theta_prev = None
            for iDC,DC in enumerate(detected_rays):
                launch_theta = DC.getLaunchVector().getTheta()/radiopropa.rad
                if iDC == (len(detected_rays)-1) or iDC == 0:
                    if iDC == 0: 
                        launch_lower.append(launch_theta-step)
                    if iDC == (len(detected_rays)-1): 
                        launch_upper.append(launch_theta+step)
                elif abs(launch_theta - launch_theta_prev) > 1.1*step: ##take 1.1 times the step to be sure the next ray is not in the bundle of the previous one
                    launch_upper.append(launch_theta_prev+step)
                    launch_lower.append(launch_theta-step)
                else:
                    pass
                launch_theta_prev = launch_theta

            roots = []
            #we minimize the cotangens of the zenith to reflect the same resolution in z to the different angles (vertical vs horizontal) 
            for i in range(len(launch_lower)):
                roots.append(optimize.brentq(delta_z, a=cot(launch_lower[i]), b=cot(launch_upper[i]), xtol=self.__ztol))
                if roots[-1].success :
                    theta = arccot(roots[-1].x)
                    detected_theta.append(theta)
                    detected_rays.append(shoot_ray(theta))

            t2 = time.time()

        self._rays = detected_rays
        self._results = [{'reflection':0,'reflection_case':1} for ray in detected_rays]
        self.__used_method = 'hybrid minimizer'

  
        elif self._config['propagation']['radiopropa']['mode'] == 'hybrid minimizing':
            has_reflec = (hasattr(self._medium,'reflection') or self.__ice_model_nuradio.reflection is not None)
            if isinstance(self._medium, medium_base.IceModelSimple) and not has_reflec:
                self.raytracer_hybrid_minimizer(n_reflections=self._n_reflections)
                results = []
                for iS in range(len(self._rays)):
                    results.append({'type':self.get_solution_type(iS), 
                                    'reflection':0,
                                    'reflection_case':1})
                self._results = results
            else:
                self.__logger.error("at the moment the RadioPropa ray tracer in minimize mode can only handle simple ice models without a reflective bottom")


