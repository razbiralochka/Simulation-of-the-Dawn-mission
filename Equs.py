import numpy as np
import scipy
import math as m
class calcs_class():
        def __init__(self, X,Y):
            self.CalcX = list()
            self.CalcY = list()
            self.angle_list = list()
            self.acc_list = list()
            self.t_angle = m.atan2(Y,X)+2*np.pi
            self.t_radius = m.sqrt(X * X + Y * Y)
            self.t_X = X
            self.t_Y = Y
            self.rVt = -2771.4710101908663
            self.phiVt = 26365.12487425152
            print("g_phi ", self.t_angle)
        def calc(self, args):
            Pr = args[0]
            Pp = args[3]
            Pvr = args[1]
            Pvp = args[2]
            Pm = -1
            SunParam = 132712440018 * pow(10, 9);
            self.CalcX.clear()
            self.CalcY.clear()
            self.angle_list.clear()
            self.acc_list.clear()
            F0 = 3 * 90 * 0.001
            m0 = 747.1 + 425 + 45.6
            a0 = F0/m0
            c0 = 26000
            X = 1.494819803027423E+08
            Y = 1.166705563891368E+07
            radius = m.sqrt(X * X + Y * Y) * 1000
            scale = radius

            angle = m.atan2(Y,X)

            V0 = np.array([-3.267334027266680E+00, 3.345804245220296E+01]) * 1000
            mV0 = m.sqrt(V0[0] * V0[0] + V0[1] * V0[1])

            rV = (X * V0[0] + Y * V0[1]) / radius
            phiV = m.sqrt(mV0 * mV0 - rV * rV)

            radius /= scale

            rV /= np.sqrt(SunParam/scale)
            phiV /= np.sqrt(SunParam / scale)
            c0 /= np.sqrt(SunParam / scale)
            T = 510*24*3600/scale/np.sqrt(scale/SunParam)

            a0 = a0 / (SunParam/scale**2)
            a = a0
            mas = 0
            t = 0

            dt = T/10000
            while t < T:
                self.CalcX.append(radius * m.cos(angle))
                self.CalcY.append(radius * m.sin(angle))

                foo = m.atan2(Pvr, Pvp)

                self.angle_list.append(m.degrees(foo))


                s = Pvp / np.sqrt(Pvp ** 2 + Pvr ** 2)
                c = Pvr / np.sqrt(Pvp ** 2 + Pvr ** 2)

                flag = (Pm/c0+np.sqrt(Pvr**2+Pvp**2)/(1-mas)) > 0


                self.acc_list.append(flag)
                dr = rV
                dphi = (phiV/radius)
                dVr = (pow(phiV, 2)/radius-1/pow(radius, 2)+c*flag*a0/(1-mas)*(1/radius)**2)
                dVphi = (-(rV*phiV)/radius+s*flag*a0/(1-mas)*(1/radius)**2)
                dm = (a0 * flag) / (c0*radius**2)

                dPr =Pp*phiV/pow(radius,2) + Pvr*(pow(phiV,2)/pow(radius,2)-2/pow(radius,3))\
                     -Pvp*rV*phiV/radius**2+2*a0*flag/pow(radius,3)*(Pm/c0+np.sqrt(Pvr**2+Pvp**2)/(1-mas))
                dPp = 0
                dPvr = -Pr+Pvp*(phiV/radius)
                dPvp = (Pvp*rV-2*Pvr*phiV-Pp)/radius
                dPm = -a0*flag/pow(radius,2)*np.sqrt(Pvr**2+Pvp**2)/(1-mas)**2

                radius += dr*dt
                angle += dphi*dt
                rV += dVr*dt
                phiV += dVphi*dt
                mas += dm*dt
                H = dr * Pr + dphi * Pp + dVr * Pvr + dVphi * Pvp + dm * Pm
                Pr += dPr*dt
                Pp += dPp*dt
                Pvr += dPvr*dt
                Pvp += dPvp*dt
                Pm += dPm*dt

                #print(H)
                t += dt
            self.rVt /= np.sqrt(SunParam/scale)
            self.phiVt /= np.sqrt(SunParam/scale)
            err = np.array([radius-self.t_radius, angle-self.t_angle, Pvr, Pvp])
            return err

        def get_points(self):
            return self.CalcX, self.CalcY, self.angle_list, self.acc_list


        def fit(self ,args):
            inp = args
            jac = np.zeros((4,4))
            err = 1000
            i = 1
            while err > 0.05:
                k = 1
                d = 1e-0
                h = d * np.diag(np.full(4, 1))
                i+=1
                out = -self.calc(inp)
                prev_err = np.linalg.norm(self.calc(inp))
                print("iter: ", i, " err ", prev_err)
                print(out)
                print(inp)
                jac[:, 0] = (self.calc(inp + h[0]) - self.calc(inp)- h[0]) / 2 * d
                jac[:, 1] = (self.calc(inp + h[1]) - self.calc(inp)- h[1]) / 2 * d
                jac[:, 2] = (self.calc(inp + h[2]) - self.calc(inp)- h[2]) / 2 * d
                jac[:, 3] = (self.calc(inp + h[3]) - self.calc(inp)- h[3]) / 2 * d



                LU = scipy.linalg.lu_factor(jac)
                darg = scipy.linalg.lu_solve(LU, out)
                inp += k*darg/np.linalg.norm(darg)
                err = np.linalg.norm(self.calc(inp))
                if err > prev_err:
                    inp -= k * darg / np.linalg.norm(darg)
                    k/=2
                    inp += k * darg / np.linalg.norm(darg)

                if k < 1e-5:
                    d /= 2
                    break

            print(inp)
            return inp