import Gadget as G
import numpy as np
import pylab as plt

#Units
Gravitational_constant = 6.6740831e-11 # m**3/(kg.s**2)
Kilo_parsec            = 3.0857e19     # m
Solar_mass             = 1.98855e30    # kg

#Transformation factor
factor = Kilo_parsec/Solar_mass

#Simu data
snap = '032'
flav = 'BH2015'
agn_good = [0,0,0,1,1,1,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,1]
direct_progenitor_032 = [0,0,0,0,36,32,82,0,0,0,0,0,2,0,0,0,183,0,0,0,0,0,0,0,0,0,0,0,0,1]
direct_progenitor_041 = [0,0,0,0,0,54,160,231,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]

#Vectors declaration
r500 = np.zeros(10)
m500 = np.zeros(10)
r200 = np.zeros(10)
m200 = np.zeros(10)

 
#Calculating Rho_critic
def Snap_header(_path):
    header            = G.snapshot_header(snap_path)
    redshift_         = header.redshift
    hubble_parameter_ = header.hubble
    omega_matter_     = header.omega_m
    omega_lambda_     = header.omega_l
    return redshift_, hubble_parameter_, omega_matter_, omega_lambda_

def Particle_data(_particle_type,_path,_redshift,_hubble_parameter):
    pos_  = G.read_block(_path,'POS ', _particle_type)/((1+_redshift)*_hubble_parameter)
    mass_ = G.read_block(_path,'MASS', _particle_type)*1e10/_hubble_parameter
    return pos_, mass_

def SubHalo_data(_path):    
    fsub_ = G.read_block(_path, 'FSUB')
    grnr_ = G.read_block(_path, 'GRNR')
    #ncon_ = G.read_block(_path, 'NCON') parece que no esta definido en el pipeline de python
    spos_ = G.read_block(_path, 'SPOS')
    return fsub_, grnr_, spos_    #fsub_, grnr_, ncon_, spos_

def Rho_critic(_redshift,_hubble_parameter,_omega_matter,_omega_lambda):
    E = omega_matter*(1+_redshift)**3+omega_lambda
    rho_critic_ = factor*3*100**2*hubble_parameter**2*E/(8*np.pi*Gravitational_constant)
    print 'RHO CRITIC: ',rho_critic_
    return rho_critic_

def Distance_respect_to_r0(_pos1,_pos0,_projected):
    if _projected==False:
        distance_ = np.sqrt((_pos1[:,0]-_pos0[0])**2 + (_pos1[:,1]-_pos0[1])**2 + (_pos1[:,2]-_pos0[2])**2)
    if _projected==True:
        distance_ = np.sqrt((_pos1[:,0]-_pos0[0])**2 + (_pos1[:,1]-_pos0[1])**2)
    return distance_

def Sort_distances_and_masses_in_terms_of_distance(_distances, _masses):
    sorted_distances_, sorted_masses_ = zip(*sorted(zip(_distances,_masses)))
    return sorted_distances_, sorted_masses_

def Interpolation(_x1,_x2,_y1,_y2,_xcrit):
    a = (_y2-_y1)/(_x2-_x1)
    b = _y1-_x1*a
    ycrit_ = a*_xcrit+b
    return ycrit_


def Total_mass_inside_Rcrit(_sorted_distances, _sorted_masses, _rho_critic, _redshift, _hubble_parameter):
    i=0
    rsphere = 2.5
    msphere = 0.0
    deltar  = 2.5
    rho_sphere_array=[0]
    rho_sphere = 0
    rho_crit500=500*2346.5392#_rho_critic
    rho_crit200=200*_rho_critic

    r = _sorted_distances
    m = _sorted_masses
    flag500=False
    while True:
        if r[i]<rsphere:
            msphere=msphere+m[i]
            i+=1
        else:
            V=(4*np.pi*rsphere**3)/3.
            rho_sphere=msphere/V
            rho_sphere_array=np.append(rho_sphere_array, rho_sphere)
            if rho_sphere > rho_crit200:
                if rho_sphere < rho_crit500 and flag500==False:
                    flag500 = True
                    r500_ = Interpolation(rho_sphere_array[-2],rho_sphere_array[-1],rsphere-deltar,rsphere,rho_crit500)
                    m500_ = ((4/3.)*np.pi*r500_**3)*rho_crit500
                rsphere=rsphere+deltar
            else:
                r200_ = Interpolation(rho_sphere_array[-2],rho_sphere_array[-1],rsphere-deltar,rsphere,rho_crit200)
                m200_ = ((4/3.)*np.pi*r200_**3)*rho_crit200
                break
        if i>len(r):
            print 'NO MORE CHANCES'
            break
    return r500_, m500_, r200_, m200_

def BCG_masses(_pos_s, _mass_s, _origin, _r500, _redshift,_hubble_parameter):
    stars_distances_respect_origin = Distance_respect_to_r0(_pos_s, _origin, False)
    sorted_stars_distances, sorted_stars_masses = Sort_distances_and_masses_in_terms_of_distance(stars_distances_respect_origin, _mass_s)

    m30_     = 0
    m50_     = 0
    m70_     = 0
    m10R500_ = 0

    r30     = 30
    r50     = 50
    r70     = 70
    r10R500 = 0.1*_r500

    for i in range(len(_mass_s)):
        if sorted_stars_distances[i]<r30:
            m30_ = m30_ + sorted_stars_masses[i]
        if sorted_stars_distances[i]<r50:
            m50_ = m50_ + sorted_stars_masses[i]
        if sorted_stars_distances[i]<r70:
            m70_ = m70_ + sorted_stars_masses[i]
        if sorted_stars_distances[i]<r10R500:
            m10R500_ = m10R500_ + sorted_stars_masses[i]
        if sorted_stars_distances[i] > np.max([r30,r50,r70 ,r10R500]):
            break
    return m30_, m50_, m70_, m10R500_


#Main Code
i=0
while i<4:
    i+=1
    if agn_good[i]:
        print '==============================================================================================================='
        region = 'D'+str(i)
        region+': I AM A GOOD REGION'
        path      = '/gss/gss_work/DRES_murante/CLUSTERS/Dianoga/D'+str(i)+'/'+flav+'/'
        snap_path = path+'snap_'+snap
        sub_path  = path+'Subfind/groups_'+snap+'/sub_'+snap
        print 'SNAP PATH ',snap_path
        print 'SUB  PATH ',sub_path
        print '-----------------------'+region+': RESULTS--------------------------'

        redshift, hubble_parameter, omega_matter, omega_lambda = Snap_header(snap_path)    #reading the snaps header
        rho_critic = Rho_critic(redshift,hubble_parameter,omega_matter,omega_lambda)       #setting cosmological parameters
        
        pos_g  , mass_g  = Particle_data(0,snap_path,redshift,hubble_parameter)            #reading positions and mass of gas         particles
        pos_dm , mass_dm = Particle_data(1,snap_path,redshift,hubble_parameter)            #reading positions and mass of dark matter particles                                                   
        pos_s  , mass_s  = Particle_data(4,snap_path,redshift,hubble_parameter)            #reading positions and mass of star        particles                                                       
        pos_bh , mass_bh = Particle_data(5,snap_path,redshift,hubble_parameter)            #reading positions and mass of black holes particles                                                                

        fsub, grnr, spos = SubHalo_data(sub_path)                                          #invoking subfind blocks

#-------loop to find the principal group in the snap       
        for j in range(10): 
            origin = spos[fsub[j],:]/((1+redshift)*hubble_parameter)     
            g_particle_distance  = Distance_respect_to_r0(pos_g  , origin, False)
            dm_particle_distance = Distance_respect_to_r0(pos_dm , origin, False)
            s_particle_distance  = Distance_respect_to_r0(pos_s  , origin, False)
            bh_particle_distance = Distance_respect_to_r0(pos_bh , origin, False)

            distances_total_particles = np.concatenate((g_particle_distance, dm_particle_distance, s_particle_distance, bh_particle_distance),axis=0)
            masses_total_particles    = np.concatenate((mass_g, mass_dm, mass_s, mass_bh),axis=0)

            sorted_distances, sorted_masses = Sort_distances_and_masses_in_terms_of_distance(distances_total_particles , masses_total_particles)
            r500[j], m500[j], r200[j], m200[j] = Total_mass_inside_Rcrit(sorted_distances, sorted_masses, rho_critic, redshift, hubble_parameter) 
            print ' '
            print '..................DATA FOR GROUP: '+str(j)+'........................'
            print ('R500 = {r5:5.5g}  M500 = {m5:10.5g}'.format(r5=r500[j],m5=m500[j]))
            print ('R200 = {r2:5.5g}  M200 = {m2:10.5g}'.format(r2=r200[j],m2=m200[j]))

 
#-------Biased Progenitor-BCGS masses
        print ' '
        print '******BIASED PROGENITOR******'

        max_index = np.where(m500==m500.max())
        main_group_array = max_index[0]
        main_group_scalar = main_group_array[0]
        r500max_array  = r500[max_index]
        r500max_scalar = r500max_array[0]
        m500max_array  = m500[max_index]
        m500max_scalar = m500max_array[0]
        r200max_array  = r200[max_index]
        r200max_scalar = r200max_array[0]
        m200max_array  = m200[max_index]
        m200max_scalar = m200max_array[0]

        m30, m50, m70, m10R500 = BCG_masses(pos_s, mass_s, origin, r500max_scalar, redshift, hubble_parameter)

        print ' '
        print 'MAIN GROUP: ', main_group_scalar
        print ('R500 = {r5:5.5g}  M500 = {m5:10.5g}'.format(r5=r500max_scalar,m5=m500max_scalar))
        print ('R200 = {r2:5.5g}  M200 = {m2:10.5g}'.format(r2=r200max_scalar,m2=m200max_scalar))
        print ('M30  = {m30:10.5g}    M50 = {m50:10.5g}    M70 = {m70:10.5g}    M10%R500 = {m10r5:10.5g}'.format(m30=m30,m50=m50,m70=m70,m10r5=m10R500))
        print ' '
        
        massesfile = open(snap+'_'+flav+'_BiasedProgenitor_Masses.txt','a')
        massesfile.write('{reg:4s} {r5:6.5g} {m5:10.5g} {r2:6.5g} {m2:10.5g} {m30:10.5g} {m50:10.5g} {m70:10.5g} {m10r5:10.5g} \n'
				 .format(reg=region,r5=r500max_scalar,m5=m500max_scalar,r2=r200max_scalar,m2=m200max_scalar,m30=m30,m50=m50,m70=m70,m10r5=m10R500))
        massesfile.close()


#-------Direct Progenitor-BCGS masses
        print ' '
        print '******DIRECT PROGENITOR******'
        print ' '

        if snap=='032':
            progenitor = direct_progenitor_032[i]
        if snap=='042': 
            progenitor = direct_progenitor_041[i]

        group         = grnr[progenitor]
        main_of_group = fsub[group]
        
        if group==main_group_scalar:
            if  main_of_group==progenitor:
                print 'DIRECT PROGENITOR IS THE SAME AS THE BIASED PROGENITOR'
                r500_progenitor     = r500max_scalar
                m500_progenitor     = m500max_scalar
                r200_progenitor     = r200max_scalar
                m200_progenitor     = m200max_scalar
                m30_progenitor      = m30
                m50_progenitor      = m50
                m70_progenitor      = m70
                m10R500_progenitor  = m10R500
            else: 
                print 'DIRECT PROGENITOR IS THE MAIN GROUP BUT IS NOT THE BIASED PROGENITOR'
                origin = spos[progenitor,:]/((1+redshift)*hubble_parameter)
                m30_progenitor, m50_progenitor, m70_progenitor, m10R500_progenitor = BCG_masses(pos_s, mass_s, origin, progenitor, redshift, hubble_parameter) 
        else:
            print 'DIRECT PROGENITOR.S GROUP IS DIFFERENT OF 0'
            r500_progenitor = r500[group]
            m500_progenitor = m500[group]
            r200_progenitor = r200[group]
            m200_progenitor = m200[group]
            origin = spos[progenitor,:]/((1+redshift)*hubble_parameter)
            m30_progenitor, m50_progenitor, m70_progenitor, m10R500_progenitor = BCG_masses(pos_s, mass_s, origin, r500_progenitor, redshift, hubble_parameter)
     
        print ' '
        print 'PROGENITOR GROUP: ', group
        print ('R500 = {r5:5.5g}  M500 = {m5:10.5g}'.format(r5=r500_progenitor,m5=m500_progenitor))
        print ('R200 = {r2:5.5g}  M200 = {m2:10.5g}'.format(r2=r200_progenitor,m2=m200_progenitor))
        print ('M30  = {m30:10.5g}    M50 = {m50:10.5g}    M70 = {m70:10.5g}    M10%R500 = {m10r5:10.5g}'.format(m30=m30_progenitor,m50=m50_progenitor,m70=m70_progenitor,m10r5=m10R500_progenitor))
        print ' '

        massesfile = open(snap+'_'+flav+'_DirectProgenitor_Masses.txt','a')
        massesfile.write('{reg:4s} {prog:6d} {grnr:6d} {fsub:6d} {r5:9.5g} {m5:10.5g} {r2:6.5g} {m2:10.5g} {m30:10.5g} {m50:10.5g} {m70:10.5g} {m10r5:10.5g} \n'
                                 .format(reg=region,prog=progenitor,grnr=group,fsub=main_of_group,r5=r500_progenitor,m5=m500_progenitor,r2=r200_progenitor,m2=m200_progenitor,m30=m30_progenitor,m50=m50_progenitor,m70=m70_progenitor,m10r5=m10R500_progenitor))
        massesfile.close()







