import numpy as np
from scipy.special import psi as digamma
import lal


#Dephasing (Eqn 5.12 - 5.20 of arXiv: 2212.13095) -----------
def psiTH_new(f,mass_1,mass_2,spin1z,spin2z):

    # Polygamma function with n=0 (digamma) -------
    def B2(x,y):
        z = x+1j*y
        temp = digamma(z)
        return temp.imag

    mass_1, mass_2 = mass_1*lal.MTSUN_SI, mass_2*lal.MTSUN_SI
    eta = (mass_1*mass_2)/(mass_1+mass_2)**2
    delta = (mass_1-mass_2)/(mass_1+mass_2)
    chi_s = (spin1z+spin2z)/2.
    chi_a = (spin1z-spin2z)/2.
    k1 = np.sqrt(1.-spin1z**2)
    k2 = np.sqrt(1.-spin2z**2)
    k_s = (k1+k2)/2.
    k_a = (k1-k2)/2.
    B21 = B2(3,2*spin1z/k1)
    B22 = B2(3,2*spin2z/k2)
    B2_s = (B21+B22)/2.
    B2_a = (B21-B22)/2.
    v = (np.pi*(mass_1+mass_2)*f)**(1./3.)
    con = 3./(128.*eta)/v**5
        
    delta_psi_5 = -(10./9.)*((1.- 3.*eta)*chi_s*(1 + 9.*chi_a**2 + 3.*chi_s**2) + delta*(1.-eta)*chi_a*(1 + 3.*chi_a**2 + 9.*chi_s**2))
    
    delta_psi_5_l = 3.*delta_psi_5
    
    delta_psi_7 = (5./168.)*(delta*chi_a*(-1667. - 4371.*chi_a**2 - 13113.*chi_s**2 + 616.*eta**2*(1 + 3.*chi_a**2 + 9.*chi_s**2) + 5.*eta*(311. + 807.*chi_a**2 +2421.*chi_s**2)) + chi_s*(840.*eta**2*(9.*chi_a**2 + 3.*chi_s**2 + 1.) + eta*(38331.*chi_a**2 + 12777.*chi_s**2 + 4889.) - 13113.*chi_a**2 - 4371.*chi_s**2 - 1667.))
    
    delta_psi_8_a = -(5./27.)*(144.*np.pi*delta*(eta-1.)*chi_a**3 + 48.*np.pi*delta*(eta-1)*chi_a + 3.*(278.*eta**2 - 370.*eta + 75.)*chi_a**4 + (-36.*eta**2 + 213.*eta - 67.)*chi_a**2 + chi_s*(-12.*delta*(eta**2 + 190.*eta - 75.)*chi_a**3 + 2.*delta*(10.*eta**2 + 124.*eta - 67.)*chi_a + 432.*np.pi*(3.*eta -1)*chi_a**2 + 48.*np.pi*(3.*eta - 1)) + chi_s**2*(432.*np.pi*delta*(eta -1)*chi_a + 90.*(36.*eta**2 - 62.*eta + 15.)*chi_a**2 - 172.*eta**2 + 303.*eta - 67.) + chi_s**3*(12.*delta*(21.*eta**2 - 130.*eta + 75.)*chi_a + 144.*np.pi*(3.*eta - 1)) + 3.*(82.*eta**2 -250.*eta + 75.)*chi_s**4 - 12.*(2.*eta**2 - 4.*eta + 1.))
    
    delta_psi_8_b = -(20./9.) * ((delta*(2.*eta - 1.)*k_a + (-1. - 2.*eta**2 + 4.*eta)*k_s)*(1. + 6.*chi_a**4 + 13.*chi_s**2 + 6.*chi_s**4 + chi_a**2*(13. + 36.*chi_s**2)) - 2.*(k_a*(1. - 4.*eta + 2.*eta**2) + delta*(1. - 2.*eta)*k_s)*chi_a*chi_s*(13. + 12.*(chi_a**2 + chi_s**2)))
    
    delta_psi_8_c = (80./9.)*(B2_s*((2.*eta**2 - 4.*eta + 1.)*chi_s*(9.*chi_a**2 + 3.*chi_s**2 +1.) - delta*(2.*eta - 1.)*chi_a*(3.*chi_a**2 + 9.*chi_s**2 + 1.)) + B2_a*(3.*(2.*eta**2 - 4.*eta + 1.)*chi_a**3 + 9.*delta*(1. - 2.*eta)*chi_a**2*chi_s + (2.*eta**2 - 4.*eta + 1.)*chi_a*(9.*chi_s**2 + 1.) - delta*(2.*eta - 1.)*chi_s*(3.*chi_s**2 + 1.)))
    
    delta_psi_8 = delta_psi_8_a + delta_psi_8_b + delta_psi_8_c
    
    delta_psi_8_l = -3.*delta_psi_8
    
    delta_psi = con*(delta_psi_5*v**5 + delta_psi_5_l*v**5*np.log(v) + delta_psi_7*v**7 + delta_psi_8*v**8 + delta_psi_8_l*v**8*np.log(v))
    
    
    return delta_psi


#-------------- ISCO for KBH -------------- (arXiv: 2108.05861)


def f_isco_Msolar_KBH(mass_1,mass_2,spin1zz,spin2zz):
    
    def r_hat_isco(chi):
        z1 = 1 + (1 - chi**2)**(1./3.)*((1 + chi)**(1./3.) + (1 - chi)**(1./3.))
        z2 = np.sqrt(3*chi**2 + z1**2)
        if chi == 0: return 3 + z2
        else: return 3 + z2 - (chi/np.abs(chi))*np.sqrt((3 - z1)*(3 + z1 + 2*z2))
    def E_hat_isco(chi):
        return np.sqrt(1 - (2.)/(3*r_hat_isco(chi)))
    def L_hat_isco(chi):
        return (2./(3*np.sqrt(3.)))*(1 + 2*np.sqrt(3*r_hat_isco(chi) - 2))

    k01 = -1.2019
    k02 = -1.20764
    k10 = 3.79245
    k11 = 1.18385
    k12 = 4.90494
    zeta = 0.41616
    k00 = -3.821158961
    mass_1 = mass_1*lal.MTSUN_SI
    mass_2 = mass_2*lal.MTSUN_SI
    M = mass_1 + mass_2 # initial total mass
    eta = (mass_1*mass_2)/(mass_1 + mass_2)**2

    S_hat = (spin1zz*mass_1**2 + spin2zz*mass_2**2)/(M**2*(1 - 2*eta))
    Erad_by_M = (0.0559745*eta + 0.580951*eta**2 - 0.960673*eta**3 + 3.35241*eta**4)\
             *((1 + S_hat*(-0.00303023 - 2.00661*eta + 7.70506*eta**2))/(1 + S_hat*(-0.067144 - 1.47569*eta + 7.30468*eta**2)))
    M_f = M*(1 - Erad_by_M)  # mass of the final black hole
    a_tot = (spin1zz*mass_1**2 + spin2zz*mass_2**2)/(mass_1 + mass_2)**2
    a_eff = a_tot + zeta*eta*(spin1zz + spin2zz)
    chi_f = a_tot + eta*(L_hat_isco(a_eff) - 2*a_tot*(E_hat_isco(a_eff) - 1))\
          + (k00 + k01*a_eff + k02*a_eff**2)*eta**2 + (k10 + k11*a_eff + k12*a_eff**2)*eta**3  # final spin of the KBH
    
    omega_hat_isco = 1./((r_hat_isco(chi_f))**(3./2.) + chi_f)
    
    return omega_hat_isco/(np.pi*M_f)


## Source model for bnary compact objects with arbitrary tidal heating 
####################################################################################
def binary_compact_object(frequency_array,mass_1,mass_2,spin1z,spin2z,luminosity_distance,
                          theta_jn,dH,**kwargs):
    """_summary_

    Args:
        frequency_array (array_like): frequency array
        mass_1 (float): mass of the first compact object in solar masses
        mass_2 (float): mass of the second compact object in solar masses
        spin1z (float): dimensionless spin of the first compact object
        spin2z (float): dimensionless spin of the second compact object
        luminosity_distance (float): luminosity distance in Mpc
        theta_jn (float): inclination angle in radians
        dH (float): tidal heating parameter [H = (1 + dH)]
                    (0 for black holes, -1 for neutron stars/perfectly reflecting compact objects)
        **waveform_kwargs: additional keyword arguments for the waveform model: minimum_frequency, maximum_frequency.
    
    Returns:
        dict: dictionary containing the waveform data
    """
    maximum_frequency = f_isco_Msolar_KBH(mass_1, mass_2, spin1z, spin2z)
    waveform_kwargs = dict(minimum_frequency=20.0, maximum_frequency=maximum_frequency)
    waveform_kwargs.update(kwargs)

    frequency_bounds = ((frequency_array >= waveform_kwargs['minimum_frequency']) *
                        (frequency_array <= waveform_kwargs['maximum_frequency']))
    
    hp, hc = binary_compact_object_waveform(frequency_array=frequency_array, 
                                            mass_1=mass_1, mass_2=mass_2, 
                                            spin1z=spin1z, spin2z=spin2z, luminosity_distance=luminosity_distance, 
                                            theta_jn=theta_jn, 
                                            dH=dH)
    hp *= frequency_bounds
    hc *= frequency_bounds

    return {"plus": hp, "cross": hc}


## TaylorF2_TidalHeated with 4.5PN NS phase (arXiv: 2304.11185) + 
# 3.5PN spinning phase + tidal heating contribution (arXiv: 2212.13095)
########################################################################

def binary_compact_object_waveform(frequency_array,
                        mass_1=10., mass_2=10., 
                        spin1z=0., spin2z=0., 
                        luminosity_distance=200., theta_jn=0.,
                        dH=0.):
    """_summary_

    Args:
        f (array_like): frequency array
        mass_1 (float, optional): Mass of the first compact object in solar masses. Defaults to 10..
        mass_2 (float, optional): Mass of the second compact object in solar masses. Defaults to 10..
        spin1z (float, optional): Dimensionless spin of the first compact object. Defaults to 0..
        spin2z (float, optional): Dimensionless spin of the second compact object. Defaults to 0..
        luminosity_distance (float, optional): Luminosity distance in Mpc. Defaults to 200..
        tc (float, optional): Coalescence time in seconds. Defaults to 0..
        phic (float, optional): Coalescence phase in radians. Defaults to 0..
        theta_jn (float, optional): Inclination angle in radians. Defaults to 0..
        dH (float, optional): Tidal heating parameter [H = (1 + dH)]
                    (0 for black holes, -1 for neutron stars/perfectly reflecting compact objects)
                    Defaults to 0..
    Returns:
        array_like, array_like: array of plus and cross polarisation strain modes 
                                evaluated at the input frequency array
    """

    PI = np.pi
    log = np.log
    sin = np.sin
    cos = np.cos
    GammaE = 0.577215664901532
    
    # convert to sec
    Ms = mass_1+mass_2
    M = Ms * lal.MTSUN_SI
    DLt = luminosity_distance * 1e6 * lal.PC_SI / lal.C_SI

    # get sym and asym chi combinations
    chi_s = 0.5*(spin1z+spin2z)
    chi_a = 0.5*(spin1z-spin2z)

    '''
    Mc is in sec, e.g., Mc = 10*MTSUN_SI (for 10 solar mass)
    DL is in sec, e.g., DL = 100*1e6*PC_SI/C_SI (for 100 Mpc)
    '''
    frequency_array[frequency_array == 0] = 1e-20 # to avoid divide by zero
    
    eta = (mass_1*mass_2)/(mass_1+mass_2)**2
    delta = (1.-4.*eta)**0.5
    v  = (PI*M*frequency_array)**(1./3.)
    # v[v == 0] = 1e-20 # to avoid error in log(v) terms
    
    beta = (113./12.)*(chi_s + delta*chi_a - (76.*eta/113.)*chi_s)    
    sigma = chi_a**2*((81./16.) - 20*eta) + (81.*chi_a*chi_s*delta)/8. + chi_s**2*((81./16.) - eta/4.)    
    epsilon = (502429./16128. - 907.*eta/192.)*delta*chi_a + (5.*eta**2/48. - 73921.*eta/2016. + 502429./16128.)*chi_s
    gamma = (732985./2268. - 24260.*eta/81. - 340.*eta**2/9.)*chi_s + (732985./2268. + 140.*eta/9.)*delta*chi_a
    
    
    # Amplitude function (without TH contribution)
    ###################################################
    # Leading order amplitude ------------------   
    def A0(mass_1,mass_2,DLt):
        eta = mass_1*mass_2/(mass_1+mass_2)**2
        M = (mass_1+mass_2)*lal.MTSUN_SI
        A = ((5.*eta/24.)**0.5/np.pi**(2./3.))*(M**(5./6.)/DLt)
        return A

    A = A0(mass_1,mass_2,DLt)
    
    a0 = 1.
    
    a1 = 0.
    
    a2 = 11.*eta/8. + 743./672.
    
    a3 = beta/2. - 2.*PI
    
    a4 = 1379.*eta**2/1152. + 18913.*eta/16128. + 7266251./8128512. - sigma/2.
    
    a5 = 57.*PI*eta/16. - 4757.*PI/1344. + epsilon
    
    a6 = 856.*GammaE/105. + 67999.*eta**3/82944. - 1041557.*eta**2/238048. - 451.*PI**2*eta/96. + 10.*PI**2/3. + 3526813753.*eta/27869184. - 29342493702821./500716339200. + 856.*np.log(4*v)/105.
    
    a7 = -1349.*PI*eta**2/24192. - 72221.*PI*eta/24192. - 5111593.*PI/2709504.

    #---------------------------------------------------------------------------------------------
    amp = A*frequency_array**(-7./6.)*(a0 + v*a1 + v**2*a2 + v**3*a3 + v**4*a4 + v**5*a5 + v**6*a6 + v**7*a7)
    #---------------------------------------------------------------------------------------------
    
    # Phase function (with TH contribution)
    ###################################################
    # 3.5PN phasing (point particle limit)
    p0 = 1.

    p1 = 0.

    p2 = (3715./756. + (55.*eta)/9.)

    p3 = (-16.*PI + (113.*delta*chi_a)/3. + (113./3. - (76.*eta)/3.)*chi_s)

    p4 = (15293365./508032. + (27145.*eta)/504.+ (3085.*eta**2)/72. + (-405./8. + 200.*eta)*chi_a**2 - (405.*delta*chi_a*chi_s)/4. + (-405./8. + (5.*eta)/2.)*chi_s**2)


    p5 = (38645.*PI/756. - 65.*PI*eta/9. - gamma)

    p5L = p5*3*log(v)

    p6 = (11583231236531./4694215680. - 640./3.*PI**2 - 6848./21.*GammaE + eta*(-15737765635./3048192. + 2255./12.*PI**2) + eta*eta*76055./1728. - eta*eta*eta*127825./1296. \
         - (6848./21.)*log(4.) + PI*(2270.*delta*chi_a/3. + (2270./3. - 520.*eta)*chi_s) + (75515./144. - 8225.*eta/18.)*delta*chi_a*chi_s \
         + (75515./288. - 263245.*eta/252. - 480.*eta**2)*chi_a**2 + (75515./288. - 232415.*eta/504. + 1255.*eta**2/9.)*chi_s**2)

    p6L = -(6848./21.)*log(v)

    p7 = (((77096675.*PI)/254016. + (378515.*PI*eta)/1512.- (74045.*PI*eta**2)/756. + (-25150083775./3048192. + (10566655595.*eta)/762048. - (1042165.*eta**2)/3024. + (5345.*eta**3)/36.
         + (14585./8. - 7270.*eta + 80.*eta**2)*chi_a**2)*chi_s + (14585./24. - (475.*eta)/6. + (100.*eta**2)/3.)*chi_s**3 + delta*((-25150083775./3048192.
         + (26804935.*eta)/6048. - (1985.*eta**2)/48.)*chi_a + (14585./24. - 2380.*eta)*chi_a**3 + (14585./8. - (215.*eta)/2.)*chi_a*chi_s**2)))

    # 4PN + 4.5PN nonspinning contribution (arXiv:2304.11185)

    p8L = (-2550713843998885153./276808510218240. + 90490.*PI**2/189. + 36812.*GammaE/63. + 1011020.*log(2.)/1323. + 78975.*log(3.)/196. + 18406.*log(v)/63. \
            +(680712846248317./42247941120. - 109295.*PI**2/224. + 3911888.*GammaE/1323. + 9964112.*log(2.)/1323. - 78975.*log(3.)/49. + 1955944.*log(v)/1323.)*eta \
            +(-7510073635./3048192. + 11275.*PI**2/144.)*eta**2 - 1292395.*eta**3/12096. + 5975.*eta**4/96.)*log(v)

    p9 = (105344279473163./18776862720. - 640.*PI**2/3. - 13696.*GammaE/21. - 13696.*log(4*v)/21. \
            + (-1492917260735./134120448. + 2255.*PI**2/6.)*eta + 45293335.*eta**2/127008. + 10323755.*eta**3/199584.)*PI

    #-----------------------------------------------------------------------------------------------------
    phase = (3./(128.*v**5*eta))*(p0 + v*p1 + v**2*p2 + v**3*p3+ v**4*p4 + v**5*(p5+p5L) + v**6*(p6+p6L) + v**7*p7 + v**8*p8L + v**9*p9)
    #-----------------------------------------------------------------------------------------------------
    # Add tidal heating contribution to the phase
    phase += (1. + dH) * psiTH_new(frequency_array,mass_1,mass_2,spin1z,spin2z)
    #-----------------------------------------------------------------------------------------------------   
    
    hp = 0.5*(1+(cos(theta_jn))**2)*amp*(cos(phase) - 1j*sin(phase))
    hc = -1j*cos(theta_jn)*amp*(cos(phase) - 1j*sin(phase))

    return hp, hc
    
