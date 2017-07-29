!***********************************************************************
!
subroutine calcscm(fi,gradfi,ifi)
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use coef
  use coefb
  use variables
  use buoy
  use time_mod
  use gradients
  use omega_turb_models
  use fieldmanipulation ! volume weighted average function

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(prec), dimension(nxyza) :: fi
  real(prec), dimension(3,nxyza) :: gradfi

  !
  ! local variables
  !
  integer ::    i, j, k, inp
  real(prec) :: gam, prtr, apotime, const, urfrs, urfms, &
                utp, vtp, wtp, utn, vtn, wtn, &
                dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                genp, genn, sut, &
                uttbuoy, vttbuoy, wttbuoy, &
                etarlzb,etarng,reta
  real(prec) :: ssq, sqrt2r
  real(prec) :: f_1, f_2, f_mu,re_y, wdis! za low-re k-epsilon modele
  real(prec) :: tmp,alphast,domegap,domegan,vist ! za k-omega modele
  real(prec) :: lengthrke, lengthles, lengthdes, dimmax, desdissiptke ! des model

  sqrt2r = 1./dsqrt(2.0d0)

  ! variable specific coefficients:
  idir=ifi
  gam=gds(ifi)
  prtr=prtinv(ifi)

  ! calculate gradient: 
  if (lstsq) then
    call grad_lsq_qr(fi,gradfi,2)
  elseif (gauss) then
    call grad_gauss(fi,gradfi(1,:),gradfi(2,:),gradfi(3,:))
  endif

  !
  ! calculate source terms integrated over volume
  if(ifi.eq.ite) then

  ! calculate production ... gen(ij)
  call find_strain_rate
  !
  do k=3,nkmm
  do i=3,nimm
  do j=3,njmm

  inp=lk(k)+li(i)+j

  if(levm) then
  !=========================================================
  ! standard production
  !=========================================================
  ssq=strain(inp)*strain(inp)
  gen(inp)=abs(vis(inp)-viscos)*ssq


  !=====production limiter for sst and sas models:=======================
  ! 10*bettainf=10*0.09=0.9 -> see below todo bettast za low-re
  if (sst.or.sas) then
  ! high-re version.....................................................
  gen(inp)=min(gen(inp),0.9*den(inp)*te(inp)*ed(inp))                !
  if (lowre) then
  ! low-re version of wilcox and sst k-omega.............................
  tmp=10.*0.09*(4./15.+(den(inp)*te(inp)/(8.*viscos*ed(inp)))**4)  &  !
              /(1.    +(den(inp)*te(inp)/(8.*viscos*ed(inp)))**4)     !  
  gen(inp)=min(gen(inp),tmp*den(inp)*te(inp)*ed(inp))                 !
  ! ....................................................................!
  end if 
  end if

  if(rng) then
  !=======================================================================
  ! strain = sqrt (sij*sij) za rng k-epsilon ! fixme da li treba da se racuna ge sa ovim strain!????
  ! s = sqrt (2*sij*sij) za realizable k-epsilon !!!
  !=======================================================================
   strain(inp) = strain(inp) * sqrt2r
  !=======================================================================
  endif

  !
  else if(lasm) then
  !=========================================================
  ! exact production (iskoristi gradijente koje vec imas dudx = gradu(1,inp)),...
  !=========================================================

  dudx = gradu(1,inp)
  dudy = gradu(2,inp)
  dudz = gradu(3,inp)

  dvdx = gradv(1,inp)
  dvdy = gradv(2,inp)
  dvdz = gradv(3,inp)

  dwdx = gradw(1,inp)
  dwdy = gradw(2,inp)
  dwdz = gradw(3,inp)

  gen(inp)=-den(inp)*(uu(inp)*dudx+uv(inp)*(dudy+dvdx)+ &
                      uw(inp)*(dudz+dwdx)+vv(inp)*dvdy+ &
                      vw(inp)*(dvdz+dwdy)+ww(inp)*dwdz)

  end if !![production calculation ]

  enddo
  enddo
  enddo

  !
  !----------------------------
  !  call calcstress
  !  call calcheatflux
  !----------------------------

  do k=3,nkmm
  do i=3,nimm
  do j=3,njmm

  inp=lk(k)+li(i)+j
  !
  genp=max(gen(inp),zero)
  genn=min(gen(inp),zero)
  !
  !=====================================
  ! [source terms]: isothermal
  !=====================================
  ! add production term to the rhs (same for all models):
  su(inp)=genp*vol(inp)              

  ! add destruction term to the lhs:

  if(stdkeps.or.durbin.or.rng.or.realizable) then
  !======================================================================
  sp(inp)=ed(inp)*den(inp)*vol(inp)/(te(inp)+small)

  !  detached eddy simulation
  if(ldes) then 
    ! k-eps length scale
    lengthrke = te(inp)**1.5_dp/ed(inp)

    ! dmax=max(dx,dy,dz)
    dimmax = max ( (x(inp)-x(inp-nj)), (y(inp)-y(inp-1)), (z(inp)-z(inp-nij)) ) 

    ! les length scale, cdes = 0.61
    lengthles = 0.61_dp*dimmax 

    ! des lengthscale
    lengthdes = min(lengthrke,lengthles) 

    ! des dissipation of tke
    desdissiptke = den(inp)*te(inp)**1.5/lengthdes

    ! add negative rhs source term to diagonal on lhs systrem matrix and divide by the unknown.
    sp(inp) = desdissiptke*vol(inp)/(te(inp)+small)
  endif   

  ! move negative prodction to lhs diagonal divided by the unknown.
  sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)  
  !======================================================================
  endif

  if(wilcox.or.sst.or.sas.or.earsm_wj.or.earsm_m)then
  !======================================================================
  ! [source terms]: isothermal
  ! note there is possibility to add a source term to eliminate 
  ! non-physical decay of turbulence variables in the freestream
  ! for external aerodynamic problems
  ! reference:
  ! spalart, p. r. and rumsey, c. l., "effective inflow conditions for 
  ! turbulence models in aerodynamic calculations," aiaa journal,
  ! vol. 45, no. 10, 2007, pp. 2544 - 2553.
  !======================================================================
  ! add sustain terms (only for sst!):
  !  su(inp)=su(inp)+0.09*tein*edin*den(inp)*vol(inp)
  ! high-re version.....................................................
  sp(inp)=bettast*ed(inp)*den(inp)*vol(inp)    
  if(lowre) then                                        
  ! low-re version of wilcox and sst k-omega.............................
    tmp =   0.09*(4./15.+(den(inp)*te(inp)/(8.*viscos*ed(inp)))**4) & !
                /(1.    +(den(inp)*te(inp)/(8.*viscos*ed(inp)))**4)   !
    sp(inp)=tmp*ed(inp)*den(inp)*vol(inp)                             !       
  ! ....................................................................!
  endif
  ! if production terms are negative, move them to lhs:
  sp(inp)=sp(inp)-genn*vol(inp)/te(inp)
  !======================================================================
  endif
  ! end: add destruction term to the lhs

  !
  !=====================================
  ! unsteady term
  !=====================================
  if(bdf) then
  !=======================================================================
  !    three level implicit time integration method:
  !    in case that btime=0. --> implicit euler
  !=======================================================================
  apotime=den(inp)*vol(inp)/timestep
  sut=apotime*((1+btime)*teo(inp)-0.5*btime*teoo(inp))
  su(inp)=su(inp)+sut
  sp(inp)=sp(inp)+apotime*(1+0.5*btime)
  !=======================================================================
  endif
  !
  !=====================================
  ! [source terms]: buoyancy
  !=====================================
  if(lcal(ien).and.lbuoy) then
  !----
  if(boussinesq) then
     uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)*beta
     vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)*beta
     wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)*beta
  else !if(boussinesq.eq.0) then
     uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)/(t(inp)+273.)
     vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)/(t(inp)+273.)
     wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)/(t(inp)+273.)
  end if
  !----
  utp=max(uttbuoy,zero)
  vtp=max(vttbuoy,zero)
  wtp=max(wttbuoy,zero)
  utn=min(uttbuoy,zero)
  vtn=min(vttbuoy,zero)
  wtn=min(wttbuoy,zero)
  !----
  su(inp)=su(inp)+utp+vtp+wtp
  sp(inp)=sp(inp)-(utn+vtn+wtn)/(te(inp)+small)
  !----
  end if

  ! end of ite volume source terms
  enddo
  enddo
  enddo

  !****************************************
  elseif(ifi.eq.ied) then
  !****************************************

  !  magstrain =  volumeweightedaverage(graded)
  !  write(66,*)'volume weighted average of dissipation grad.',magstrain

  !======================================================================
  ! terms for menter sst model:
  if (sst.or.sas) call calculate_stuff_for_sst
  ! terms for earsm_wj and earsm_m model:
  if (earsm_wj.or.earsm_m) call calculate_stuff_for_earsm
  !======================================================================

  do k=3,nkmm
  do i=3,nimm
  do j=3,njmm

  inp=lk(k)+li(i)+j
  !
  genp=max(gen(inp),zero)
  genn=min(gen(inp),zero)
  !
  !===============================================
  ! [hybrid alpha model: ce2=ce1+0.48/alpha
  !===============================================
  if (alphamodel.and.durbin) then
  al_rans(inp)=te(inp)**(3./2.)/(ed(inp)+small)
  al_les(inp)=vol(inp)**(1./3.)
  alph(inp)=max(1.d0,al_rans(inp)/al_les(inp))
  c2=c1+0.48/alph(inp)
  endif

  !======================================================================
  ! [hybrid alpha model: ce2=ce1+0.48/alpha, samo k-omega sst verzija!
  ! experimentalna varijanta!
  ! ovde je diff je funkcija polozaja i ekvivalent je 0.48 kod k-eps
  !  0.48=c2-c1 => u sst modelu: diff(inp)=bettasst(inp)-alphasst(inp)
  ! aktivira se tako sto sst=true and sas=true.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(sst.and.alphamodel) then
  al_rans(inp)=te(inp)**(1/2.)/(ed(inp)*cmu)
  al_les(inp)=vol(inp)**(1/3.)
  !al_kolm(inp)=((viscos**3.)/(ed(inp)*te(inp)*cmu))**(1/4.)
  alph(inp)=max(1.d0,al_rans(inp)/al_les(inp))
  !alph(inp)=max(1.d0,0.5*al_rans(inp)/al_les(inp))
  bettasst(inp)=alphasst(inp)+diff(inp)/alph(inp)
  !======================================================================
  endif


  !=====cross diffusion for sst model====================================
  !if (sst.or.sas.or.earsm_wj.or.earsm_m) then
    domegap=max(domega(inp),zero)
    domegan=min(domega(inp),zero)
  !end if
  !======================================================================

  !
  !=====================================
  ! [source terms]: isothermal
  !=====================================
  sv(inp)=0.0d0 ! this is why we use sv in bcscalarm.f!!! so there we add only gent calculated in wallbc.f

  !
  !=====standard k-epsilon or keps+durbin=============================
  if (stdkeps.or.durbin) then
  su(inp)=c1*genp*ed(inp)*vol(inp)/(te(inp)+small)
  sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
  sp(inp)=sp(inp)-c1*genn*vol(inp)/(te(inp)+small)
  !=====end: standard k-epsilon or keps+durbin=============================
  endif 

  !=====realizable k-epsilon==============================================
  if (realizable) then

  genp=max(strain(inp),zero)
  genn=min(strain(inp),zero)

  etarlzb = strain(inp)*te(inp)/(ed(inp)+small)
  c1 = max(0.43,etarlzb/(etarlzb+5.))
  su(inp)=c1*genp*ed(inp)*vol(inp)
  sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/ &
                           ( te(inp)+sqrt(viscos/densit*ed(inp)) )
  sp(inp) = sp(inp) - c1*genn*ed(inp)*vol(inp)
  !=====end:realizable k-epsilon==========================================
  endif
  !
  !=====renormalization group (rng) k-epsilon=============================
  if (rng) then
  su(inp)=c1*genp*ed(inp)*vol(inp)/(te(inp)+small)
  sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
  sp(inp)=sp(inp)-c1*genn*vol(inp)/(te(inp)+small)

  etarng = strain(inp)*te(inp)/(ed(inp)+small)
  reta = cmu*etarng**3*(1-etarng/4.38)/(1+0.012*etarng**3)
  if (etarng.le.4.38d0) then
  sp(inp)=sp(inp)+reta*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
  else
  su(inp)=su(inp)-reta*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
  end if
  !=====end:renormalization group (rng) k-epsilon=========================
  endif

  !
  !=====low-re lam-bremhorst k-epsilon====================================
  if (lowre_lb) then
  !
  ! lam-bremhorst :
  !

  wdis = walldistance(inp)
  re_y = den(inp)*sqrt(te(inp))*wdis/(viscos+small)
  ret(inp)=den(inp)*te(inp)**2/(viscos*ed(inp)+small) 
  f_mu = (1.-exp(-0.0165*re_y))**2*(1.+20.5/ret(inp))
  f_1 = 1.+(0.05/f_mu)**3
  f_2 = 1.-exp(-ret(inp)**2) 
  ! 
  ! launder-sharma :
  !
  !ret(inp)=den(inp)*te(inp)**2/(viscos*ed(inp)+small) 
  !f_1 = 1.0d0
  !f_2 = 1.0d0-0.3d0*exp(-ret(inp)**2)
  !
  ! for all:
  su(inp) = f_1*c1*genp*ed(inp)*vol(inp)/(te(inp)+small)
  sp(inp) = f_2*c2*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
  sp(inp) = sp(inp)-f_1*c1*genn*vol(inp)/(te(inp)+small)
  ! za low-re modele - premesti f_mu u ret - da znas u modvis kada vidis ret to je zapravo f_mu  
  ret(inp) = f_mu
  !=====end:low-re lam-bremhorst k-epsilon================================
  endif ! low-re lam-bremhorst  


  !
  if (wilcox) then
  !=====wilcox  k-omega===================================================
  su(inp)=alpha*genp*ed(inp)*vol(inp)/te(inp)
  if (lowre) then
  ! low-re version......................................................
  ! let's find alpha*                                                  !
  alphast=(0.024+(densit*te(inp))/(6.*viscos*ed(inp)))  &            !
         /(1.+(densit*te(inp))/(6.*viscos*ed(inp)))                  !
  tmp=alpha/alphast*(1./9.+(densit*te(inp))/(2.95*viscos*ed(inp))) & !
                    /((1.+(densit*te(inp))/(2.95*viscos*ed(inp))))   !
  su(inp)=tmp*genp*ed(inp)*vol(inp)/te(inp)                          !
  ! ...................................................................!
  endif

  ! add destruction term to the lhs:
  sp(inp)=betta*den(inp)*ed(inp)*vol(inp) 
  sp(inp)=sp(inp)-alpha*genn*vol(inp)/(te(inp)+small)
  !=====end: wilcox  k-omega===============================================
  end if

  !=====menter k-omega sst and sas ========================================
  ! ovde ce sad uzeti limitiran gen(inp) jer je gore usao u if sst
  ! ovde \alpha nije konstanta kao u std vec se racuna u "calc_stuff_for_sst",
  ! takodje ovde deli produkciju sa \nu_t umesto da mnozi sa \omega/k.
  if (sst) then
  vist = (vis(inp)-viscos)/densit
  su(inp)=alphasst(inp)*genp*vol(inp)/vist
  su(inp)=su(inp)+domegap*vol(inp)
  ! add sustain terms
  !  su(inp)=su(inp)+bettasst(inp)*edin*edin*den(inp)*vol(inp)
  ! add destruction term to the lhs:
  vist = (vis(inp)-viscos)/densit
  sp(inp)=bettasst(inp)*den(inp)*ed(inp)*vol(inp) 
  sp(inp)=sp(inp)-alphasst(inp)*genn*vol(inp)  &
                   /(vist*ed(inp))
  !  sp(inp)=sp(inp)-domegan*vol(inp)/ed(inp)
  !=====end: menter k-omega sst ============================================
  end if

  ! sas only. adds qsas production
  if (sas) then
  !=====sas only ===============================================
  vist = (vis(inp)-viscos)/den(inp)
  su(inp)=alphasst(inp)*genp*vol(inp)/vist
  su(inp)=su(inp)+domegap*vol(inp)
  ! add sustain term
  !  su(inp)=su(inp)+bettasst(inp)*edin*edin*den(inp)*vol(inp)
  ! add sas production term (!) :
  su(inp)=su(inp)+qsas(inp)*vol(inp)
  ! add destruction term to the lhs:
  vist = (vis(inp)-viscos)/densit
  sp(inp)=bettasst(inp)*den(inp)*ed(inp)*vol(inp) 
  sp(inp)=sp(inp)-alphasst(inp)*genn*vol(inp)  &
                   /(vist*ed(inp))
  !  sp(inp)=sp(inp)-domegan*vol(inp)/ed(inp)
  !=====end: sas only ===============================================
  end if


  if (earsm_wj.or.earsm_m) then
  !=====earsm_wj and earsm_m models standalone or in sas hybrid version ==========
  su(inp)=alphasst(inp)*genp*ed(inp)*vol(inp)/(te(inp)+small)   !earsm
  su(inp)=su(inp)+domegap*vol(inp)
  ! add sustain terms
  !  su(inp)=su(inp)+bettasst(inp)*edin*edin*den(inp)*vol(inp)

  ! sas also adds qsas production
  if (sas) then
  ! add sas production term
  su(inp)=su(inp)+qsas(inp)*vol(inp)
  end if

  ! add destruction term to the lhs:
  sp(inp)=bettasst(inp)*den(inp)*ed(inp)*vol(inp) 
  sp(inp)=sp(inp)-alphasst(inp)*genn*vol(inp) &
                   /(te(inp)+small)              !earsm
  sp(inp)=sp(inp)-domegan*vol(inp)/ed(inp)
  !=====end: earsm_wj and earsm_m models standalone or in sas hybrid version =====
  end if

  !
  !=====================================
  ! unsteady term
  !=====================================
  if(bdf) then
  !=======================================================================
  !    three level implicit time integration method:
  !    in case that btime=0. --> implicit euler
  !=======================================================================
  apotime=den(inp)*vol(inp)/timestep
  sut=apotime*((1+btime)*edo(inp)-0.5*btime*edoo(inp))
  su(inp)=su(inp)+sut
  sp(inp)=sp(inp)+apotime*(1+0.5*btime)
  !=======================================================================
  endif
  !
  !=====================================
  ! [source terms]: buoyancy
  !=====================================
  if(lcal(ien).and.lbuoy) then
  const=c3*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
  !----
  if(boussinesq) then
     uttbuoy=-gravx*utt(inp)*const*beta
     vttbuoy=-gravy*vtt(inp)*const*beta
     wttbuoy=-gravz*wtt(inp)*const*beta
  else !if(boussinesq.eq.0) then
     uttbuoy=-gravx*utt(inp)*const/(t(inp)+273.)
     vttbuoy=-gravy*vtt(inp)*const/(t(inp)+273.)
     wttbuoy=-gravz*wtt(inp)*const/(t(inp)+273.)
  end if
  !----
  utp=max(uttbuoy,zero)
  vtp=max(vttbuoy,zero)
  wtp=max(wttbuoy,zero)
  utn=min(uttbuoy,zero)
  vtn=min(vttbuoy,zero)
  wtn=min(wttbuoy,zero)
  !----
  su(inp)=su(inp)+utp+vtp+wtp
  sp(inp)=sp(inp)-(utn+vtn+wtn)/(ed(inp)+small)
  !----
  end if

  ! end of ied volume source terms
  enddo
  enddo
  enddo
  !--------------------------------------
  end if

  !
  ! calculate terms integrated over surfaces
  ! only inner surfaces
  ! east cell - face
  call fluxscm(nj,1,nij,fi,gradfi,ifi, &
               ar1x,ar1y,ar1z, &
               fx,ae,aw,f1)

  ! north cell - face
  call fluxscm(1,nij,nj,fi,gradfi,ifi, &
               ar2x,ar2y,ar2z, &
               fy,an,as,f2)

  ! top   cell - face
  call fluxscm(nij,nj,1,fi,gradfi,ifi, &
               ar3x,ar3y,ar3z, &
               fz,at,ab,f3)

  urfrs=urfr(ifi)
  urfms=urfm(ifi)

  do k=3,nkmm
  do i=3,nimm
  do j=3,njmm

  inp=lk(k)+li(i)+j

  if(cn) then
  ! crank-nicolson stuff - only once:
  ae(inp)=0.5d0*ae(inp)
  an(inp)=0.5d0*an(inp)
  at(inp)=0.5d0*at(inp)
  aw(inp)=0.5d0*aw(inp)
  as(inp)=0.5d0*as(inp)
  ab(inp)=0.5d0*ab(inp)
  ! crank-nicolson time stepping source terms
  apotime=den(inp)*vol(inp)/timestep
  if(ifi.eq.ite) then
  su(inp)=su(inp)+(ae(inp)*teo(inp+nj)  + aw(inp)*teo(inp-nj)+    &
                   an(inp)*teo(inp+1)   + as(inp)*teo(inp-1)+     &
                   at(inp)*teo(inp+nij) + ab(inp)*teo(inp-nij))+  &
          (apotime-ae(inp)-aw(inp)                                &
                  -an(inp)-as(inp)                                &
                  -at(inp)-ab(inp))*teo(inp) 
  else ! ifi==ied   
  su(inp)=su(inp)+(ae(inp)*edo(inp+nj)  + aw(inp)*edo(inp-nj)+    &
                   an(inp)*edo(inp+1)   + as(inp)*edo(inp-1)+     &
                   at(inp)*edo(inp+nij) + ab(inp)*edo(inp-nij))+  &
          (apotime-ae(inp)-aw(inp)                                &
                  -an(inp)-as(inp)                                &
                  -at(inp)-ab(inp))*edo(inp)      
  endif ! checking if this is tke or epsilon field      
  sp(inp)=sp(inp)+apotime
  ! end of crank-nicolson time stepping source terms
  ! end of crank-nicolson stuff:
  endif

  ! main diagonal term assembly
  ap(inp)=ae(inp)+aw(inp)+an(inp)+as(inp)+at(inp)+ab(inp)+sp(inp)
  ! underelaxation:
  ap(inp)=ap(inp)*urfrs
  su(inp)=su(inp)+urfms*ap(inp)*fi(inp)
  enddo
  enddo
  enddo
  !
  ! solving linear system:
  call sipsol(fi,ifi)
  !call cgstab_sip(fi,ifi)

  ! these field values cannot be negative
  if(ifi.eq.ite.or.ifi.eq.ied) then
     do inp=icst,icen
     fi(inp)=max(fi(inp),small)
     enddo
  endif

return
end subroutine
