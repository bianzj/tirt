"""Opt-in finite-cone ports for whole-crown optical transport.

No foliage-volume cubature. Ground-radius and angular quadrature average
geometric routes from whole crowns. The routes condition already escaped
foliage flux; they are not raw solid-angle fractions or fitted BRF gains.
Other mountain interactions retain the existing terrain view factors.
"""
import math
from functools import lru_cache
import numpy as np

from base import optical_crown as crown
from base.optical_crown_profile import source_profile, escape_ratios, beam_down_escape
from base.gap import gap_probability_crown_analytical as gap
from base.optical_terrain import exchange_geometry
from base.terrain_rt import _shape_values, _projected_gap
from base.solar import slope_direction
from base import optical_crown_intrinsic as intrinsic


def _intrinsic(geo):
    return len(geo)>4 and geo[4]=='intrinsic_crown'


def _profile(shape,rho,tau,order,geo):
    return (intrinsic.source_profile(shape,rho,tau,order) if _intrinsic(geo) else
            source_profile(shape,rho,tau,order))


def _primary(shape,solar,rho,tau,soil,order,mode,scales,geo):
    return (intrinsic.primary(shape,solar,rho,tau,soil,order,mode,scales) if _intrinsic(geo) else
            crown.primary(shape,solar,rho,tau,soil,order,mode,'beam_diffusion_profile',scales))


def cone_hit(origins, direction, height, radius):
    """Nearest forward finite-cone side hit, excluding its buried base."""
    k = radius/height
    a = direction[0]**2+direction[1]**2-k*k*direction[2]**2
    b = 2*(origins[:,:2]@direction[:2]+k*k*(height-origins[:,2])*direction[2])
    c = (origins[:,:2]**2).sum(axis=1)-k*k*(height-origins[:,2])**2
    discriminant = b*b-4*a*c
    answer = np.full(len(origins),np.inf)
    if abs(a)<1e-12:
        roots = [-c/np.where(abs(b)>1e-12,b,np.nan)]
    else:
        roots = [(-b-sign*np.sqrt(np.maximum(discriminant,0)))/(2*a) for sign in (-1,1)]
    for t in roots:
        z = origins[:,2]+t*direction[2]
        valid = (discriminant>=0)&(t>1e-7)&(z>=0)&(z<=height)
        answer = np.minimum(answer,np.where(valid,t,np.inf))
    return answer


@lru_cache(maxsize=65536)
def route(height,radius,centre,crown_radius,mu,azimuth):
    """Fractions to own slope soil, outside flat forest, and open sky.

    mu is signed relative to the slope normal. Crowns are uniformly sampled
    per horizontal area, as in the original scene. The radius quadrature is
    over ground positions, never over positions inside a crown.
    """
    nodes,weights = np.polynomial.legendre.leggauss(16)
    radii = radius*np.sqrt((nodes+1)/2)
    weights = weights/2
    beta = math.atan2(height,radius)
    sn,cs = math.sin(beta),math.cos(beta)
    normal = np.array([sn,0.,cs])
    tangent = np.array([cs,0.,-sn])
    direction = mu*normal+math.sqrt(max(0.,1-mu*mu))*(
        math.cos(azimuth)*tangent+math.sin(azimuth)*np.array([0.,1.,0.]))
    origins = np.c_[radii+centre*sn,np.zeros(len(radii)),
                    height*(1-radii/radius)+centre*cs]
    # Own-crown attenuation is already represented by the analytical escape
    # kernel. Route only after leaving the whole-crown bounding sphere.
    origins += crown_radius*direction
    side = cone_hit(origins,direction,height,radius)
    flat_ground = np.full(len(radii),np.inf)
    if direction[2]<-1e-12:
        flat_ground = -origins[:,2]/direction[2]
    terrain = np.minimum(side,flat_ground)
    bottom,top = centre-crown_radius,centre+crown_radius
    # First entry to the outside flat-canopy envelope. Its response includes
    # its own gaps and soil, so no additional foliage interception multiplier.
    entry = np.full(len(radii),np.inf)
    candidates = [np.zeros(len(radii))]
    if abs(direction[2])>1e-12:
        candidates.extend([(level-origins[:,2])/direction[2] for level in (bottom,top)])
    a = float(direction[:2]@direction[:2])
    if a>1e-12:
        b = 2*(origins[:,:2]@direction[:2])
        c = (origins[:,:2]**2).sum(axis=1)-radius*radius
        disc = b*b-4*a*c
        candidates.append(np.where(disc>=0,(-b+np.sqrt(np.maximum(disc,0)))/(2*a),np.inf))
    for t in candidates:
        point = origins+np.where(np.isfinite(t),t,0)[:,None]*direction
        valid = ((t>=-1e-8)&np.isfinite(t)&(point[:,2]>=bottom-1e-7)&
                 (point[:,2]<=top+1e-7)&((point[:,:2]**2).sum(axis=1)>=radius*radius-1e-6))
        entry = np.minimum(entry,np.where(valid,np.maximum(t,0),np.inf))
    outside = np.isfinite(entry)&(entry<terrain)
    sky = ~outside&~np.isfinite(terrain)
    own = ~outside&~sky
    return tuple(float(weights@label) for label in (own,outside,sky))


@lru_cache(maxsize=2048)
def ports(shape,height,radius,centre,crown_radius,attenuation,order,depth_mode='effective_forest'):
    depth = (intrinsic.intrinsic_depth(shape) if depth_mode=='intrinsic_crown' else
             source_profile(shape,0.,0.,order)[0])
    mu,w,az = crown.quadrature(order,48)
    sums = np.zeros(3)
    total = 0.
    for m,weight in zip(mu,w):
        view = math.degrees(math.acos(m))
        kernel = 2*m*weight*(1-float(gap(*shape[:4],view)))/shape[0]
        kernel *= escape_ratios(view,depth,attenuation)[1]
        total += kernel
        for angle in az:
            sums += kernel/len(az)*np.asarray(route(height,radius,centre,crown_radius,
                                                   -float(m),math.radians(float(angle))))
    return tuple(sums/total)


@lru_cache(maxsize=4096)
def state(shape,solar,rho,tau,soil,order,mode,scales,geo):
    first,escaped,secondary,scale = _primary(shape,solar,rho,tau,soil,order,mode,scales,geo)
    matrix,exits = (intrinsic.exchange_matrix(shape,rho,tau,soil,order) if _intrinsic(geo) else
                   crown.exchange_matrix(shape,rho,tau,soil,order,'beam_diffusion_profile'))
    depth,attenuation,eup,edn = _profile(shape,rho,tau,order,geo)
    depth_mode = 'intrinsic_crown' if _intrinsic(geo) else 'effective_forest'
    rg,rf,rs = ports(shape,*geo[:4],attenuation,order,depth_mode)
    rg1,rf1,rs1 = ports(shape,*geo[:4],.5/math.cos(math.radians(solar)),order,depth_mode)
    d1 = first[1]*(intrinsic.beam_down_escape(shape,solar,order) if _intrinsic(geo) else
                  beam_down_escape(shape,solar,order))
    matrix = matrix.copy(); secondary = secondary.copy()
    matrix[2,1] *= rg
    secondary[2] *= rg1
    higher = np.linalg.solve(np.eye(3)-matrix,secondary)
    down = higher[1]*edn
    up = float(escaped.sum()+eup*higher[0]+exits[2]*higher[2])
    outside = float(d1*rf1+down*rf)
    sky = float(d1*rs1+down*rs)
    escape = np.array([eup,edn*(rf+rs),exits[2]])
    absorbed = (1-first.sum()+first.sum()-escaped.sum()-d1*(rf1+rs1)-secondary.sum()+
                (1-matrix.sum(axis=0)-escape)@higher)
    if min(up,outside,sky,absorbed)<-1e-9 or abs(up+outside+sky+absorbed-1)>1e-9:
        raise ValueError('Finite-cone port energy balance failed.')
    return higher,scale,np.array([up,outside,sky]),float(absorbed),float(d1)


def parameters(settings):
    shape = tuple(float(settings[k]) for k in ('lai','stand_density','height','width','hspot'))
    order = int(settings.get('optical_crown_quadrature',12))
    mode = settings.get('optical_crown_geometry','clumping_dual')
    if mode!='clumping_dual' or not 4<=order<=32:
        raise ValueError('Finite-cone ports require clumping_dual and quadrature in [4,32].')
    scales = (float(settings.get('optical_hotspot_inner',shape[4])),
              float(settings.get('optical_hotspot_outer',shape[3]/shape[2]))) if mode=='clumping_dual' else None
    return shape,order,mode,scales


@lru_cache(maxsize=32768)
def beam(shape,solar,view,relative,rho,tau,soil,order,mode,scales,geo):
    higher,scale,_,_,_ = state(shape,solar,rho,tau,soil,order,mode,scales,geo)
    first = crown._single(shape,solar,view,relative,rho,tau,soil,mode,scales)
    first[0] *= scale
    depth,att,_,_ = _profile(shape,rho,tau,order,geo)
    gv = float(gap(*shape[:4],view))
    return float(first.sum()+higher[0]*(1-gv)/shape[0]*escape_ratios(view,depth,att)[0]+higher[2]*gv)


@lru_cache(maxsize=4096)
def sky(shape,view,rho,tau,soil,order,mode,scales,geo):
    mu,w,az = crown.quadrature(order)
    return sum(2*m*weight/len(az)*beam(shape,float(np.degrees(np.arccos(m))),view,float(a),
               rho,tau,soil,order,mode,scales,geo) for m,weight in zip(mu,w) for a in az)


@lru_cache(maxsize=512)
def diffuse_ports(shape,rho,tau,soil,order,mode,scales,geo):
    mu,w,_ = crown.quadrature(order)
    return sum(2*m*weight*state(shape,float(np.degrees(np.arccos(m))),rho,tau,soil,
                              order,mode,scales,geo)[2] for m,weight in zip(mu,w))


@lru_cache(maxsize=32768)
def down_beam(shape,solar,view,rho,tau,soil,order,mode,scales,geo):
    higher,_,_,_,_ = state(shape,solar,rho,tau,soil,order,mode,scales,geo)
    first = _primary(shape,solar,rho,tau,soil,order,mode,scales,geo)[0]
    dep,att,_,_ = _profile(shape,rho,tau,order,geo)
    kernel = (1-float(gap(*shape[:4],view)))/shape[0]
    return kernel*(first[1]*escape_ratios(view,dep,.5/math.cos(math.radians(solar)))[1]+
                   higher[1]*escape_ratios(view,dep,att)[1])


@lru_cache(maxsize=4096)
def down_sky(shape,view,rho,tau,soil,order,mode,scales,geo):
    mu,w,_ = crown.quadrature(order)
    return sum(2*m*weight*down_beam(shape,float(np.degrees(np.arccos(m))),view,rho,tau,soil,
                                  order,mode,scales,geo) for m,weight in zip(mu,w))


def response(settings,solar,view,relative,rho,tau,soil,f=0.,sunlit=True,geo=None):
    if geo is None:
        if settings.get('optical_terrain_source_depth')=='intrinsic_crown':
            return intrinsic.reflectance(settings,solar,view,relative,rho,tau,soil,f,sunlit)
        return crown.crown_reflectance(settings,solar,view,relative,rho,tau,soil,f,sunlit)
    shape,order,mode,scales = parameters(settings)
    return ((1-f)*beam(shape,solar,view,relative,rho,tau,soil,order,mode,scales,geo) if sunlit else 0.)+(
        f*sky(shape,view,rho,tau,soil,order,mode,scales,geo) if f else 0.)


def reflectance(settings,sza,saa,vza,vaa,rho,tau,soil,diffuse_fraction=0.,sunlit=True,return_budget=False):
    """Three outgoing ports: facet-up, flat-up, facet-to-flat forest.

    Open-sky downward escape has independent crown visibility on backfaces.
    Bare mountain view factors are retained for the original upward ports.
    """
    from base.optical_terrain import forest_albedo,diffuse_albedo
    depth_mode = settings.get('optical_terrain_source_depth','effective_forest')
    if depth_mode not in ('effective_forest','intrinsic_crown'):
        raise ValueError('Unknown finite-crown source depth: '+str(depth_mode))
    if depth_mode=='intrinsic_crown':
        forest_albedo,diffuse_albedo=intrinsic.albedo,intrinsic.diffuse_albedo
    values = [sza,saa,vza,vaa,rho,tau,soil,diffuse_fraction]
    if (not np.isfinite(values).all() or min(rho,tau,soil)<0 or rho+tau>1+1e-12 or soil>1):
        raise ValueError('Invalid optical coefficients or angles.')
    shapes = _shape_values(settings)
    if len(shapes)!=1:
        raise ValueError('Experimental finite-cone ports currently require one cone class.')
    height,radius,density = shapes[0]
    if density==0:
        if depth_mode=='intrinsic_crown':
            return intrinsic.reflectance(settings,sza,vza,abs((vaa-saa+180)%360-180),
                                         rho,tau,soil,diffuse_fraction,sunlit)
        return crown.crown_reflectance(settings,sza,vza,abs((vaa-saa+180)%360-180),
                                       rho,tau,soil,diffuse_fraction,sunlit)
    if height<=0 or radius<=0 or density<0 or density*math.pi*radius**2>=1:
        raise ValueError('Invalid finite-cone geometry.')
    if not 0<=diffuse_fraction<=1 or not 0<=sza<90 or not 0<=vza<90:
        raise ValueError('Invalid illumination or direction.')
    if float(settings['lai'])<=1e-8:
        from base.optical_terrain import terrain_crown_reflectance
        return terrain_crown_reflectance(settings,sza,saa,vza,vaa,rho,tau,soil,diffuse_fraction,sunlit)
    if (min(float(settings[k]) for k in ('lai','stand_density','height','width'))<=0 or
            abs(float(settings['height'])-float(settings['width']))>1e-9):
        raise ValueError('Finite-cone ports currently require positive spherical whole crowns.')
    local = dict(settings,surface_model='plane',optical_crown_source='beam_diffusion_profile',
                 optical_crown_geometry=settings.get('optical_crown_geometry','clumping_dual'))
    # crown_bottom is normal distance from ground, not a vertical tree height.
    centre = float(settings.get('crown_bottom',settings['height']))+float(settings['height'])
    if not math.isfinite(centre) or centre<float(settings['height']):
        raise ValueError('Crown bottom must be finite and nonnegative.')
    geo = (height,radius,centre,float(settings['height']))
    if depth_mode=='intrinsic_crown':
        geo += (depth_mode,)
    beta = math.atan2(height,radius); cs = math.cos(beta)
    facet = dict(local,lai=local['lai']*cs,stand_density=local['stand_density']*cs)
    shape,order,mode,scales = parameters(facet)
    relative = abs((vaa-saa+180)%360-180)
    plane,plane_sun,_ = _projected_gap(shapes,vza,sza,relative)
    # Macro source is a hemispheric illumination, independent of the viewer.
    solar_plane = _projected_gap(shapes,sza,sza,0.)[0]
    coverage = density*math.pi*radius**2
    flat_lit = solar_plane/(1-coverage)
    factors = exchange_geometry(shapes)
    sky_flux = 1-factors.sum(axis=1)
    f = diffuse_fraction
    result = plane*((1-f)*response(local,sza,vza,relative,rho,tau,soil,0.,sunlit)*
                    (plane_sun/plane if plane else 0.)+f*sky_flux[1]*response(local,sza,vza,relative,rho,tau,soil,1.,sunlit))
    slope = math.degrees(beta)
    # The statistical cone distribution has no preferred world azimuth.
    aspects = float(saa)+np.arange(0.,360.,5.)
    cv = cs*math.cos(math.radians(vza))+math.sin(beta)*math.sin(math.radians(vza))*np.cos(np.radians(vaa-aspects))
    denom = np.maximum(cv,0).sum()
    visible_response = diffuse_response = 0.
    direct_ports = np.zeros(3)
    direct_absorbed = direct_incident = 0.
    for aspect,cosv in zip(aspects,cv):
        sv,sa,lit = slope_direction(sza,saa,slope,aspect)
        vv,va,_ = slope_direction(vza,vaa,slope,aspect)
        raa = abs((va-sa+180)%360-180)
        sv = min(sv,89.9)
        incident = max(0.,math.cos(math.radians(sv))/math.cos(math.radians(sza))) if lit and sunlit else 0.
        if incident:
            packet = state(shape,sv,rho,tau,soil,order,mode,scales,geo)
            direct_ports += incident*packet[2]/len(aspects)
            direct_absorbed += incident*packet[3]/len(aspects)
            direct_incident += incident/len(aspects)
        if cosv>1e-9:
            weight = cosv/denom
            visible_response += weight*(1-f)*incident*response(facet,sv,min(vv,89.9),raa,rho,tau,soil,geo=geo)
            sk = response(facet,sv,min(vv,89.9),raa,rho,tau,soil,1.,geo=geo)
            diffuse_response += weight*sk
        elif cosv<-1e-9:
            # Emission below a facet's tangent can still reach the global sky.
            # Legacy slope1 folds zeniths into [0,90]; recover the downward
            # hemisphere from the signed dot product, not from its return.
            downview = math.degrees(math.acos(min(1.,-float(cosv))))
            # slope1 azimuth is defined in its local frame. Route() uses its
            # outward tangent; world dot products avoid a frame convention.
            world = np.array([math.sin(math.radians(vza))*math.cos(math.radians(vaa-aspect)),
                              math.sin(math.radians(vza))*math.sin(math.radians(vaa-aspect)),
                              math.cos(math.radians(vza))])
            tangent = np.array([cs,0.,-math.sin(beta)])
            az = math.atan2(world[1],float(world@tangent))
            visible_sky = route(*geo[:4],float(cosv),az)[2]
            if visible_sky:
                direct_down = down_beam(shape,sv,downview,rho,tau,soil,order,mode,scales,geo)
                # Absolute projected area for foliage; the soil stays filtered.
                visible_response += (-cosv)/denom*(1-f)*incident*visible_sky*direct_down
                diffuse_response += (-cosv)/denom*visible_sky*down_sky(
                    shape,downview,rho,tau,soil,order,mode,scales,geo)
    up_sky,flat_sky,sky_down = diffuse_ports(shape,rho,tau,soil,order,mode,scales,geo)
    af = diffuse_albedo(local,rho,tau,soil)
    areas = np.array([coverage/cs,1-coverage])
    area_ratio = areas[0]/areas[1]
    # Column sources are Uf, U0, Df. Df is intercepted by the outside flat
    # canopy, whose full response includes all of its own scattering orders.
    operator = np.array([[up_sky*factors[0,0],up_sky*factors[0,1],0.],
                         [af*factors[1,0],0.,af*area_ratio],
                         [flat_sky*factors[0,0],flat_sky*factors[0,1],0.]])
    source = (1-f)*np.array([direct_ports[0],flat_lit*forest_albedo(local,sza,rho,tau,soil) if sunlit else 0.,direct_ports[1]])
    source += f*np.array([up_sky*sky_flux[0],af*sky_flux[1],flat_sky*sky_flux[0]])
    outgoing = np.linalg.solve(np.eye(3)-operator,source)
    incoming_facet = factors[0]@outgoing[:2]
    incoming_flat = factors[1,0]*outgoing[0]+area_ratio*outgoing[2]
    result += (1-plane)*(visible_response+diffuse_response*(f*sky_flux[0]+incoming_facet))
    result += plane*response(local,sza,vza,relative,rho,tau,soil,1.,sunlit)*incoming_flat
    if not return_budget:
        return float(result)
    direct_flat = forest_albedo(local,sza,rho,tau,soil)
    external_in = (1-f)*(areas[0]*direct_incident+areas[1]*flat_lit*bool(sunlit))+f*(areas@sky_flux)
    direct_abs = (1-f)*(areas[0]*direct_absorbed+areas[1]*flat_lit*bool(sunlit)*(1-direct_flat))
    direct_abs += f*(areas[0]*sky_flux[0]*(1-up_sky-flat_sky-sky_down)+areas[1]*sky_flux[1]*(1-af))
    absorbed = direct_abs+areas[0]*(1-up_sky-flat_sky-sky_down)*incoming_facet+areas[1]*(1-af)*incoming_flat
    open_sky = (areas[0]*outgoing[0]*sky_flux[0]+areas[1]*outgoing[1]*sky_flux[1]+
                areas[0]*((1-f)*direct_ports[2]+sky_down*(f*sky_flux[0]+incoming_facet)))
    return float(result),dict(incident=float(external_in),escaped=float(open_sky),absorbed=float(absorbed),
        residual=float(open_sky+absorbed-external_in),macro_operator=operator.tolist(),
        macro_outgoing=outgoing.tolist(),macro_incoming=[float(incoming_facet),float(incoming_flat)],
        areas=areas.tolist(),facet_direct_ports=direct_ports.tolist(),facet_diffuse_ports=[float(up_sky),float(flat_sky),float(sky_down)])
