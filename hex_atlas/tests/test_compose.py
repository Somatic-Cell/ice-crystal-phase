from pathlib import Path
import math
import sys
import unittest
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
from atlas_math import (canonical_frame, fold_direction, interpolation, gaussian_tilt_nodes,
    orientation_rule, projected_area, coupling, compose_mass, frame, make_cdf, ry, rz)


def fixture():
    it,ip,nt,np_=6,3,4,12
    area=np.zeros((it+1,ip+1));mass=np.zeros((it+1,ip+1,nt,np_))
    pattern=np.arange(1,nt*np_+1,dtype=float).reshape(nt,np_);pattern/=pattern.sum()
    for i in range(it+1):
        for j in range(ip+1):
            area[i,j]=projected_area(.1,.2,canonical_frame(i*math.pi/2/it,j*math.pi/6/ip)[:,2])
            mass[i,j]=area[i,j]*pattern
    return SimpleNamespace(it=it,ip=ip,nt=nt,np=np_,area=area,mass=mass,residual=np.zeros_like(area),radius=.1,length=.2)


class Tests(unittest.TestCase):
    def test_fold_full_sphere(self):
        for theta in np.linspace(0,math.pi,21):
            for phi in np.linspace(0,2*math.pi,43):
                k=canonical_frame(theta,phi)[:,2];t,p,q=fold_direction(k)
                self.assertTrue(0<=t<=math.pi/2 and 0<=p<=math.pi/6)
                np.testing.assert_allclose(q@canonical_frame(t,p)[:,2],k,atol=2e-14)
                np.testing.assert_allclose(q.T@q,np.eye(3),atol=2e-14)

    def test_gaussian_cdf_nodes(self):
        for sigma in (.1,3,10,90):
            n=16;b=gaussian_tilt_nodes(n,sigma);s=math.radians(sigma)
            e=math.erf(math.pi/(2*math.sqrt(2)*s))
            c=np.array([(1+math.erf((x-math.pi/2)/(math.sqrt(2)*s))/e)/2 for x in b])
            np.testing.assert_allclose(c,(np.arange(n)+.5)/n,atol=6e-14,rtol=0)
            np.testing.assert_allclose(b+b[::-1],math.pi,atol=2e-15)
        np.testing.assert_array_equal(gaussian_tilt_nodes(4,0),[math.pi/2]*4)

    def test_rule_no_rng_and_weights(self):
        config={'quadrature':{'azimuth':4,'tilt':8,'spin':3},'components':[
            {'kind':'horizontal_gaussian','sigma_deg':3,'fraction':.9},{'kind':'isotropic','fraction':.1}]}
        a,w,m=orientation_rule(config);b,v,_=orientation_rule(config)
        np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(w,v)
        self.assertFalse(m['randomness']);self.assertAlmostEqual(math.fsum(w[:96]),.9)
        self.assertAlmostEqual(math.fsum(w[96:]),.1)

    def test_bad_rule(self):
        for sigma in (-1,float('nan')):
            with self.assertRaises(ValueError):gaussian_tilt_nodes(10,sigma)
        with self.assertRaises(ValueError):orientation_rule({'quadrature':{'azimuth':2,'tilt':2,'spin':2},'components':[{'kind':'isotropic','fraction':.5}]})

    def test_interpolation_endpoints_and_partition(self):
        for i in range(7):
            for j in range(4):
                terms=interpolation(i*math.pi/12,j*math.pi/18,6,3)
                self.assertAlmostEqual(sum(w for _,_,w in terms),1)
                self.assertAlmostEqual(sum(w for ii,jj,w in terms if (ii,jj)==(i,j)),1)
        for t in np.linspace(0,math.pi/2,21):
            terms=interpolation(t,.14,6,3)
            self.assertTrue(all(w>=0 for _,_,w in terms))
            self.assertAlmostEqual(sum(w for _,_,w in terms),1)

    def test_no_incident_extrapolation(self):
        with self.assertRaises(ValueError):interpolation(-.01,.1,6,3)
        with self.assertRaises(ValueError):interpolation(.1,.6,6,3)

    def test_shift_and_reflection(self):
        a=fixture();f=canonical_frame(math.pi/6,math.pi/18);fw=frame([1,0,0])
        for reflected in (False,True):
            q=np.diag([1.,1.,-1.]) if reflected else np.eye(3)
            for shift in (0,1,2.3,-.7,11.999):
                d=2*math.pi*shift/a.np
                h=np.array([[math.cos(d),-math.sin(d),0],[math.sin(d),math.cos(d),0],[0,0,1.]])
                if reflected:h=h@np.diag([1.,-1.,1.])
                r=fw@h@f.T@q.T
                w,expected,_=coupling(a,np.array([r]),np.ones(1),[1,0,0])
                result,residual=compose_mass(a,w)
                original=a.mass[2,1]
                if reflected:original=original[:,::-1]
                lower=math.floor(shift);fraction=shift-lower
                reference=(1-fraction)*np.roll(original,lower,axis=1)+fraction*np.roll(original,lower+1,axis=1)
                np.testing.assert_allclose(result,reference,atol=3e-17,rtol=1e-12)
                self.assertAlmostEqual(result.sum(),expected);self.assertEqual(residual,0)

    def test_area_weighted_mixture(self):
        a=fixture();rot=np.array([np.eye(3),rz(math.pi/2)])
        w,area,_=coupling(a,rot,np.array([.25,.75]),[1,0,0])
        out,tail=compose_mass(a,w)
        exact=.25*projected_area(.1,.2,[1,0,0])+.75*projected_area(.1,.2,[0,-1,0])
        self.assertAlmostEqual(area,exact);self.assertAlmostEqual(out.sum()+tail,exact)

    def test_residual_not_renormalized_away(self):
        a=fixture();a.mass*=.99;a.residual=.01*a.area
        w,area,_=coupling(a,np.array([ry(.213)@rz(.7)]),np.ones(1),[1,-.4,.3])
        mass,tail=compose_mass(a,w)
        self.assertAlmostEqual(mass.sum()/area,.99);self.assertAlmostEqual(tail/area,.01)

    def test_cdf_uniform_and_zero_rows(self):
        edges=(1-np.cos(np.linspace(0,math.pi,9)))/2
        for mass in (np.repeat(np.diff(edges)[:,None]/12,12,axis=1),np.eye(8,12)):
            a,b,r=make_cdf(mass,edges)
            recovered=np.diff(a)[:,None]*np.diff(b,axis=1)
            np.testing.assert_allclose(recovered.T,mass/mass.sum(),atol=2e-15)
            self.assertAlmostEqual(r['cdf_mass'],1)
            zero=np.diff(a)==0
            if zero.any():np.testing.assert_array_equal(b[zero],np.repeat(edges[None,:],zero.sum(),axis=0))

    def test_cdf_dynamic_range_and_rejection(self):
        edges=np.linspace(0,1,9);mass=np.geomspace(1e-14,1,96).reshape(8,12)
        _,_,r=make_cdf(mass,edges);self.assertLess(r['cdf_l1_mass_error'],1e-12)
        with self.assertRaises(ValueError):make_cdf(-mass,edges)
        with self.assertRaises(ValueError):make_cdf(np.zeros_like(mass),edges)

if __name__=='__main__':unittest.main(verbosity=2)
