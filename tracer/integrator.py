import numpy as np 

class RayTracer: 

    def __init__(self, rs=1.0, r0=50.0):
        self.rs = rs
        self.r0 = r0  # camera distance from object 
        self.b_crit = (3* np.sqrt(3) / 2) * rs



    def _binet_rhs(self, u, w):
        """
        Right hand side of the Binet ODE system
        split into two ODEs of 
        du/dphi = w 
        dw/dphi = -u + (3/2) * rs * u^2 
        """

        return w, -u + 1.5 * self.rs * u**2 
    
    def _rk4_step(self, u, w, dphi):
        """
        Advance (u, w) by one RK4 step of size dphi
        u = position 
        w = velocity 
        """
        k1u, k1w = self._binet_rhs(u,w)
        k2u, k2w = self._binet_rhs(u + 0.5*dphi*k1u, w + 0.5*dphi*k1w)
        k3u, k3w = self._binet_rhs(u + 0.5*dphi*k2u, w + 0.5*dphi*k2w)
        k4u, k4w = self._binet_rhs(u + dphi*k3u, w + dphi*k3w)

        u_new = u + (dphi / 6) * (k1u + 2*k2u + 2*k3u +k4u )
        w_new = w + (dphi / 6) * (k1w + 2*k2w + 2*k3w +k4w )

        return u_new, w_new
    
    def trace_ray(self, b, dphi=0.0002, max_steps=300000):
        """
        Trace a photon with impact parameter b 

        Returns:
            captured: True if photon crossed the horizon 
            exit_phi: Total angle swept
        """

        u = 1.0 / self.r0 
        w = np.sqrt(1.0/b**2 - u**2 + self.rs * u**3)

        phi = 0.0
        passed_closet = False
        
        for _ in range(max_steps):
            if u > 1.0 / self.rs: 
                return True, phi
            
            if w < 0: 
                passed_closet = True
                
            if passed_closet and (1.0 / u) > self.r0:
                return False, phi
            
            phi += dphi
            u, w = self._rk4_step(u, w, dphi)

        return False, phi
    

    def trace_ray_path(self, b, dphi=0.0002, max_steps=300000):
        """Like trace_ray but returns the full (phi, r) trajectory."""
        u = 1.0 / self.r0
        w = np.sqrt(1.0/b**2 - u**2 + self.rs * u**3)

        phi = 0.0
        passed_closest = False
        phis, rs = [phi], [self.r0]

        for _ in range(max_steps):
            if u > 1.0 / self.rs:
                phis.append(phi)
                rs.append(self.rs)
                break

            if w < 0:
                passed_closest = True

            if passed_closest and (1.0 / u) > self.r0:
                phis.append(phi)
                rs.append(1.0 / u)
                break

            phi += dphi
            u, w = self._rk4_step(u, w, dphi)
            phis.append(phi)
            rs.append(1.0 / u)

        return np.array(phis), np.array(rs)