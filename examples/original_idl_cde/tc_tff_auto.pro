
; Routine ACCEPT_PARSE opens the ACCEPT ASCII file, extracts the profiles
;   for the cluster given by NAME, and writes those profiles to a
;   *_profiles.dat file readable by ACCEPT_READ.

pro accept_parse,name
openr,1,'all_profiles.dat'
openw,2,'ACCEPT_profiles/'+strtrim(name)+'_profiles.dat'
dummy = ' '
name1 = ' '
readf,1,name1,dummy,format="(a20,a310)"
printf,2,dummy
readf,1,name1,dummy,format="(a20,a310)"
printf,2,dummy
while not eof(1) do begin
  readf,1,name1,dummy,format="(a20,a310)"
  if name1 eq name then printf,2,dummy
  endwhile
close,1
close,2
return
end

; Routine ACCEPT_READ reads files created by ACCEPT_PARSE from the
;  master ASCII ACCEPT file.

pro accept_read,file,profiles
openr,1,file
dummy = ' '
readf,1,dummy
readf,1,dummy
nlines=1
while not eof(1) do begin
  readf,1,dummy
  nlines = nlines + 1
  endwhile
close,1
openr,1,file
readf,1,dummy
readf,1,dummy
profiles = fltarr(19,nlines-1)
readf,1,profiles
close,1
return
end

; profiles(0,*) is Rin  (Mpc)
; profiles(1,*) is Rout  (Mpc)
; profiles(2,*) is nelec  (cm^-3)
; profiles(3,*) is neerr  (cm^-3)
; profiles(4,*) is Kitpl  (keV cm^2)
; profiles(5,*) is Kflat  (keV cm^2)
; profiles(6,*) is Kerr  (keV cm^2)
; profiles(7,*) is Pitpl  (dyne cm^-2)
; profiles(8,*) is Pflat  (dyne cm^-2)
; profiles(9,*) is Perr  (dyne cm^-2)
; profiles(10,*) is Mgrav  (M_solar)
; profiles(11,*) is Merr  (M_solar)
; profiles(12,*) is Tx  (keV)
; profiles(13,*) is Txerr  (keV)
; profiles(14,*) is lambda  (ergs-cm^3/s)
; profiles(15,*) is tcool5/2  (Gyr)
; profiles(16,*) is t52err  (Gyr)
; profiles(17,*) is tcool3/2  (Gyr)
; profiles(18,*) is t32err  (Gyr)



; Routine LOGTEMP fits the logarithmic ACCEPT electron density profile
;  ln n_e (in cm^-3) to a polynomial in log r (in Mpc) of degree DEG.
;  The vector COEFFS returns the coefficients of that polynomial fit.

pro logdensity_fit,profiles,r,deg,coeffs
rin = profiles(0,*)
nbins = n_elements(rin)
rin = rin(0:nbins-1)
rout = profiles(1,*)
rout = rout(0:nbins-1)
r = (rin + rout) * 0.5
logr = alog(r)
nelec = profiles(2,*)
nelec = nelec(0:nbins-1)
lognelec = alog(nelec)
neerr = profiles(3,*)
neerr = neerr(0:nbins-1)
logneerr = alog(neerr/nelec)
yerror = logneerr
coeffs = poly_fit(logr,lognelec,deg,yerror=logneerr, $
  /double,status=ok,chisq=chi2)
print,chi2,coeffs
lognefit = 0.0 * logr
for i = 0,deg do lognefit = lognefit + coeffs(i)*logr^float(i)
nefit = exp(lognefit)
plot_oo,r,nelec,xtitle='r (Mpc)',ytitle='n!de!n (cm!u-3!n)',line=2
oplot,r,nelec+neerr,line=1
oplot,r,nelec-neerr,line=1
oplot,r,nefit,line=2
return
end


; Routine LOGTEMP fits the logarithmic ACCEPT temperature profile
;  ln kT (in keV) to a polynomial in log r (in Mpc) of degree DEG.
;  The vector COEFFS returns the coefficients of that polynomial fit.

pro logtemp_fit,profiles,r,deg,coeffs,tfit,logterr
rin = profiles(0,*)
nbins = n_elements(rin)
rin = rin(0:nbins-1)
rout = profiles(1,*)
rout = rout(0:nbins-1)
r = (rin + rout) * 0.5
logr = alog(r)
t = abs(profiles(12,*))
t = t(0:nbins-1)
logt = alog(t)
terr = profiles(13,*)
terr = terr(0:nbins-1)
logterr = alog(terr/t)
yerror = logterr
coeffs = poly_fit(logr,logt,deg,yerror=logneerr, $
  /double,status=ok,chisq=chi2)
print,chi2,coeffs
logtfit = 0.0 * logr
for i = 0,deg do logtfit = logtfit + coeffs(i)*logr^float(i)
tfit = exp(logtfit)
plot_oo,r,t,xtitle='r (Mpc)',ytitle='kT (KeV)',line=2
oplot,r,t+terr,line=1
oplot,r,t-terr,line=1
oplot,r,tfit
return
end


; Routine LOGPRESSURE fits the logarithmic ACCEPT electron pressure
;  profile ln P_e (in erg cm^-3) to a polynomial in log r (in Mpc) of
;  degree DEG.  The vector COEFFS returns the coefficients of
;  that polynomial fit.

pro logpressure_fit,profiles,r,deg,coeffs,pfit,logperr
rin = profiles(0,*)
nbins = n_elements(rin)
rin = rin(0:nbins-1)
rout = profiles(1,*)
rout = rout(0:nbins-1)
r = (rin + rout) * 0.5
logr = alog(r)
p = profiles(7,*)
p = p(0:nbins-1)
logp = alog(p)
perr = profiles(9,*)
perr = perr(0:nbins-1)
logperr = alog(perr/p)
yerror = logperr
coeffs = poly_fit(logr,logp,deg,yerror=logneerr, $
  /double,status=ok,chisq=chi2)
print,chi2,coeffs
logpfit = 0.0 * logr
for i = 0,deg do logpfit = logpfit + coeffs(i)*logr^float(i)
pfit = exp(logpfit)
plot_oo,r,p,xtitle='r (Mpc)',ytitle='P (erg cm!u-3!n)',line=2
oplot,r,p+perr,line=1
oplot,r,p-perr,line=1
oplot,r,pfit
return
end


; Routine RG_DATA determines the gravitational potential corresponding to
;  fits (PFIT and TFIT) to the ACCEPT pressure and temperature profiles.
;  RG gives the potential as the radius times the gravitational
;  acceleration in cgs units.

pro rg_data,profiles,rMpc,rg,tfit,pfit,rgerr, $
                     rMpc_fine,rg_fine,tfit_fine,pfit_fine

rin = profiles(0,*)
nbins = n_elements(rin)
rin = rin(0:nbins-1)
rout = profiles(1,*)
rout = rout(0:nbins-1)
rMpc = (rin + rout) * 0.5
ln_rMpc = alog(rMpc)
rcm = rMpc * 3.08D24
log_rMpc_fine = findgen(300)/100. - 3.
rMpc_fine = 10.0^log_rMpc_fine
ln_rMpc_fine = alog(rMpc_fine)

degt = 3
logtemp_fit,profiles,rdummy,degt,tcoeffs,tfit,logterr
kt_erg = tfit * 1.6E-9
logtfit = 0.0 * ln_rMpc_fine
for i = 0,degt do logtfit = logtfit + tcoeffs(i)*ln_rMpc_fine^float(i)
tfit_fine = exp(logtfit)
kt_erg_fine = tfit_fine * 1.6E-9
oplot,rMpc_fine,tfit_fine,line=3
stop

degp = 3
logpressure_fit,profiles,rdummy,degp,pcoeffs,pfit,logperr
logpfit = 0.0 * ln_rMpc_fine
for i = 0,degp do logpfit = logpfit + pcoeffs(i)*ln_rMpc_fine^float(i)
pfit_fine = exp(logpfit)
oplot,rMpc_fine,pfit_fine,line=3
stop

dlnp_dlnr =  0.0 * ln_rMpc
for i = 1,degp do dlnp_dlnr = dlnp_dlnr $
                          + float(i)*pcoeffs(i)*ln_rMpc^float(i-1)
dlnp_dlnr = dlnp_dlnr < (-1.0E-10)
print,dlnp_dlnr

dlnp_dlnr_fine =  0.0 * ln_rMpc_fine
for i = 1,degp do dlnp_dlnr_fine = dlnp_dlnr_fine $
                          + float(i)*pcoeffs(i)*ln_rMpc_fine^float(i-1)
dlnp_dlnr_fine = dlnp_dlnr_fine < (-1.0E-10)

mu_mp = 0.6 * 1.67D-24
rg = - kt_erg / mu_mp * dlnp_dlnr
rg_fine = - kt_erg_fine / mu_mp * dlnp_dlnr_fine
plot_oi,rMpc,rg,xtitle='r (Mpc)',ytitle = 'rg (cgs)'
oplot,rMpc_fine,rg_fine,line=3

relerr = sqrt(2.*exp(logperr)^2. + exp(logterr)^2.)
rgerr = kt_erg / mu_mp * relerr
oplot,rMpc,rg+rgerr,line=1
oplot,rMpc,rg-rgerr,line=1
stop
return
end


; Routine LAMBDA_TN03 gives the cooling function for of Tozzi &
;  Norman (2001) to the Sutherland & Dopita cooling function for
;  0.3 solar metallicity:  LAMBDA_CGS in cgs units.

pro lambda_tn03,kt_kev,lambda_cgs
c1 = 8.6d-3
c2 = 5.8d-2
c3 = 6.3d-2
alpha = -1.7
beta = 0.5
lambda = c1*kt_kev^alpha + c2*kt_kev^beta + c3
lambda_cgs = lambda * 1.0d-22
return
end

; Routine TC_TFF_AUTO gives the ratio of cooling time to free-fall time
;  derived from fits to ACCEPT profiles.

pro tc_tff_auto,name
;accept_parse,name ; THIS IS NO LONGER NEEDED!
file = 'ACCEPT_profiles/'+strtrim(name)+'_profiles.dat'
accept_read,file,profiles
rg_data,profiles,rMpc,rg,tfit,pfit,rgerr, $
                     rMpc_fine,rg_fine,tfit_fine,pfit_fine
rgmin = 2.0 * (2.5d7)^2.

rg_fine = rg_fine > rgmin
r_fine = rMpc_fine * 3.08d24
tff_fine = sqrt(2. / rg_fine) * r_fine
plot_oo,rMpc_fine,tff_fine/3.15d7,yrange=[1.0d7,1.0d11],line=1

rg = rg > rgmin
r = rMpc * 3.08d24
tff = sqrt(2. / rg) * r
oplot,rMpc,tff/3.15d7

lambda_tn03,tfit,lambda_cgs
nelec = pfit / (tfit*1.6d-9)
tc = 1.5 * (1.89 * pfit) / (nelec^2. / 1.07) / lambda_cgs
oplot,rMpc,tc/3.15d7

lambda_tn03,tfit_fine,lambda_cgs
nelec_fine = pfit_fine / (tfit_fine*1.6d-9)
tc_fine = 1.5 * (1.89 * pfit_fine) / (nelec_fine^2. / 1.07) / lambda_cgs
oplot,rMpc_fine,tc_fine/3.15d7,line=2

stop

plot_oo,rMpc_fine,tc_fine/tff_fine,yrange=[1,100],xtitle='r (Mpc)',ytitle='t!dc!n/t!dff!n',line=1
oplot,rMpc,tc/tff

trat = tc/tff
trat_fine = tc_fine/tff_fine
imin = max(where(rMpc_fine lt min(rMpc))) > 0
imax = min(where(rMpc_fine gt max(rMpc))) < n_elements(rMpc_fine)-1
minrat = min(trat_fine(imin:imax))
rmin = rMpc_fine(where(trat_fine eq minrat))

openw,1,'tc_tff_profiles/'+strtrim(name)+'_tc_tff.dat'
for i = 0,n_elements(rMpc)-1 do printf,1,rMpc(i),tc(i),tff(i),trat(i)
close,1

openw,2,'tc_tff_profiles/'+strtrim(name)+'_tc_tff_fine.dat'
for i = 0,n_elements(rMpc_fine)-1 do printf,2,rMpc_fine(i), $
               tc_fine(i),tff_fine(i),trat_fine(i)
close,2

print,name
print,"Minimum Ratio = ",minrat
print,"At Radius (in kpc) = ",rmin * 1000.0

stop

mass = rg / 6.67d-8 * r / 2.d33
plot_oo,rMpc,mass,xtitle='r (Mpc)',ytitle='M (M!dSun!n)'
oplot,rMpc,mass*(1+rgerr/rg),line=1
oplot,rMpc,mass*(1-rgerr/rg),line=1

return
end
