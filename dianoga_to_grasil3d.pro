pro dia_to_g3dv, win, ireg, auto=auto, write=write
  ;MIND THE FLAVOR, largheza box!
  ; correr con auto si se desea centrar el box usando subfind
  
  do_bh = 1
  n_s   = 0L
  n_g   = 0L
  box   = 500. ; larghezza fisica del box    
  
  simu_dir   = '/gss/gss_work/DRES_murante/CLUSTERS/Dianoga/'
  output_dir = '/pico/scratch/userexternal/cragonef/mferraro/g3d/Input/' 
 
  if ireg EQ 0 then begin
    print, 'D0 DOES NOT EXIST'
    exit
  endif
 
  regions=['none','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28','D29']

  flav='BH2015'
  ;flav='CSF2015'
  snap='091'
  
  ;ireg= 2 ; non c'e' snap_032 para el flav bhnew de la regione 1, 3, 4 

  snap_name=simu_dir+regions[ireg]+'/'+flav+'/snap_'+snap
  print, 'SNAP Name: ', snap_name
  sub_name=simu_dir+regions[ireg]+'/'+flav+'/Subfind/groups_'+snap+'/sub_'+snap
  print, 'SUB  Name: ', sub_name
  
  if keyword_set(auto)  then begin
    auto=auto-1
    readnew, sub_name, spos, 'SPOS'  ; posizioni di tutti i sottoaloni
    readnew, sub_name, gpos, 'GPOS'  ; posizioni di tutti i sottoaloni
    readnew, sub_name, smst, 'SMST'  ; masas de cada subhalo
    readnew, sub_name, fsub, 'FSUB'  ; indici del primo sottoalone di ogni alone negli array di sottoaloni (tipo spos)
    readnew, sub_name, mvir, 'MVIR'  ; masse di ogni alone
    readnew, sub_name, rvir, 'RVIR'  ; masse di ogni alone
   ; readnew, sub_name, ncon, 'NCON'  ; contaminacion de cada halo, si es distinto de cero esta contaminado 
   ; readnew, sub_name, mtot, 'MTOT'  ; contaminacion de cada halo, si es distinto de cero esta contaminado 
    
    ;  for i=0,100 do begin
    ;    masse=msub[fsub[i]:fsub[i+1]-1]
    ;    maxmas=max(masse,imax)
    ;    print,  imax
    ;  endfor
   ; print, '*********MVIR', mvir[0:10]*1e10 
   ; print, '*********MTOT', mtot[0:10]*1e10 
   ; print, '*********SPOS', spos[*,fsub[0]] 
   ; print, '*********GPOS', gpos[*,0] 
   ; print, '*********NCONT', ncon[0:10] 
    select_with_mvir = 0     ;CIN
    if select_with_mvir then begin
      mvir = mvir[where(ncon eq 0)]
      fsub = fsub[where(ncon eq 0)]
      ind=reverse(sort(mvir))
      fsub=fsub[ind] ; ora fsub e' ordinato in senso decrescente della massa dell'alone
      ;mvir=mvir[ind] ; non so se servira'
      pos_cdm=spos(*,fsub[auto])
    endif else begin         ; para la region 1 (BH2015) Mvir del cumulo principal es mas chica que la de un cumulo que esta en el borde. Subfind esta haciendo algo mal?
      mfsub = smst[4,fsub]   ; masas de los subhalos principales, aquellos que contienen BCG+ICL    ;CIN
      ind   = reverse(sort(mfsub))                                                                  ;CIN
      pos_cdm=spos(*,fsub[auto])                                                                    ;CIN
      ;pos_cdm=spos(*,89)                                                                    ;CIN
    endelse
    
    auto_str=strtrim(string(auto),2)
  endif else begin
    print, 'WARNING !!!! cmd scritto a mano'  
    pos_cdm=[492757.,498711.,497772.] ; g0016649_G/csf_w500 snap 032 passatomi da ste z=2
    auto_str=''
  endelse
  print, 'pos_cdm', pos_cdm
  
  box_str=strtrim(string(nint(box)),2)
  if do_bh then begin
    ;file_out_name=output_dir+regions[ireg]+'_'+flav+'_'+snap+'_clu'+auto_str+'_L'+box_str+'kpc_siagn.dat'    
    file_out_name=output_dir+regions[ireg]+'_'+flav+'_'+snap+'_clu'+auto_str+'_L'+box_str+'kpc_nd.dat'    
    ;file_out_name=output_dir+'prueba89.dat'    
  endif else begin
    file_out_name=output_dir+regions[ireg]+'_'+flav+'_'+snap+'_clu'+auto_str+'_L'+box_str+'kpc_noagn.dat'
  endelse  
  
  ;stop
 print, file_out_name
 
  ; test low res cioe' la ris di sempre...
  ;file_out_name="/scratch3/granato/Image/Input/dianoga_test_lr.dat"
  ;snap_name='/scratch/borgani/Dianoga/g0052436_Me14_G/csf_w500/snap_040'
  ;pos_cdm=[499631.,493374.,501427.]
  ;box_half_size=250.
  
  
  ; test low res cioe' la ris di sempre...
  ;file_out_name="/scratch3/granato/g3d/Input/g0016649.dat"
  ;snap_name='/scratch/murante/DianogaHR/g0052436_Me14_G/snap_027'
  ;pos_cdm=[499688.,493465.,501330.] ; file dianoga_test_hr.dat snap 41
  ;pos_cdm=[499971.,498070.,501836.] ; dia2_test_hr.dat snap 038
  ;pos_cdm=[499957.,498125., 497606.] ; g0016649_G/csf_bhnew-w500 snap 032 altra pos interessante: pos_cdm=[492783.,498704., 497753.]
  
  ;box_half_size=125. file dianoga_test_hr.dat snap 41 dia2_test_hr.dat snap 038
  
  
  readnew, snap_name, head, "HEAD"
  z=head.redshift
  print, 'REDSHIFT', z
  hubble=head.hubbleparam
  omega_m=head.omega0
  omega_lambda=head.omegalambda
  
  ; gas
  readnew, snap_name, cldx_g, "CLDX", parttype=0
  readnew, snap_name, pos_g, "POS ", parttype=0
  readnew, snap_name, M_g, "MASS", parttype=0
  readnew, snap_name, Z_g, "Zs  ", parttype=0
  readnew, snap_name, rho_g, "RHO ", parttype=0 ; this is rho_code:  
  readnew, snap_name, T_g, "TEMP", parttype=0
  
  ; rho_phys = rho_code * (1 + z)^3 * m_unit / l_unit^3 * hpar^2
  ; where m_unit=10^10 Msun and l_unit=kpc
  
  ;stars 
  readnew, snap_name, pos_s, "POS ", parttype=4
  readnew, snap_name, M_s, "MASS", parttype=4
  readnew, snap_name, Mi_s, "iM  ", parttype=4
  readnew, snap_name, Z_s, "Zs  ", parttype=4
  readnew, snap_name, age_s, "AGE ", parttype=4
  
  ; BHs
  if do_bh then begin
    readnew, snap_name, pos_bh, "POS ", parttype=5
    readnew, snap_name, M_bh, "BHMA", parttype=5
    ; in snapshot accretion rate given in units of 10^10 h^-1 Msun / (0.98 Gyr). 
    ; 0.98 is the ratio between 3.08 (coming from parsec in cm) and 3.155 (coming from number of seconds in one year)
    readnew, snap_name, Mdot_bh, "BHMD ", parttype=5
    ;readnew, snap_name, sl, "ACRBH", parttype=5
  endif
  
  box_half_gadget=box*(1+z)*hubble/2.
  
  ; select only particles in the desired box
  ind_box_g=where( abs(pos_g[0,*]-pos_cdm[0]) lt box_half_gadget and abs(pos_g[1,*]-pos_cdm[1]) lt box_half_gadget and abs(pos_g[2,*]-pos_cdm[2]) lt box_half_gadget)
  ind_box_s=where( abs(pos_s[0,*]-pos_cdm[0]) lt box_half_gadget and abs(pos_s[1,*]-pos_cdm[1]) lt box_half_gadget and abs(pos_s[2,*]-pos_cdm[2]) lt box_half_gadget)
  pos_g=pos_g[*,ind_box_g] & M_g=M_g[ind_box_g] & Z_g=Z_g[*,ind_box_g] & rho_g=rho_g[ind_box_g] & cldx_g=cldx_g[ind_box_g] & t_g=t_g[ind_box_g]
  pos_s=pos_s[*,ind_box_s] & M_s=M_s[ind_box_s] & Mi_s=Mi_s[ind_box_s] & Z_s=Z_s[*,ind_box_s] & age_s=age_s[ind_box_s] 
  if do_bh then begin
    ind_box_bh=where( abs(pos_bh[0,*]-pos_cdm[0]) lt box_half_gadget and abs(pos_bh[1,*]-pos_cdm[1]) lt box_half_gadget and abs(pos_bh[2,*]-pos_cdm[2]) lt box_half_gadget)
    pos_bh=pos_bh[*,ind_box_bh] & M_bh=M_bh[ind_box_bh] & Mdot_bh=Mdot_bh[ind_box_bh] 
  endif
  
  
  ; select only gas particles which can host dust. cldx_g is the fraction of mass of gas particle in cold phase
  ; if comment this to use all particles, remember also below in units conversion comment multiplications bu cldx_g
  ;ind_dusty=where(cldx_g gt 0 ) ; PRIMA VERSIONE
  ;ind_dusty=where(cldx_g gt 0 or t_g lt 1e5) ; ALLA MUPPICI..
  ind_dusty=where(cldx_g gt 0 or t_g lt 1e5, nind_dusty)
  print, 'N DUSTY= ', nind_dusty,'*********************************************************************************************' 
  ;ind_dusty=where(t_g lt 1e15) ; ALLA CAZZO
; ind_dusty=where(cldx_g gt 0 or t_g lt 1e6) ; !!!!!!!!!!!  TEST
  pos_g=pos_g[*,ind_dusty] & M_g=M_g[ind_dusty] & Z_g=Z_g[*,ind_dusty] & rho_g=rho_g[ind_dusty] & cldx_g=cldx_g[ind_dusty]
  
  for i=0,2 do begin ; refer coordinates to centre of selected object
    pos_g[i,*]=pos_g[i,*]-pos_cdm[i]
    pos_s[i,*]=pos_s[i,*]-pos_cdm[i]
    if do_bh then begin
      pos_bh[i,*]=pos_bh[i,*]-pos_cdm[i]
    endif
  endfor
  
   
  n_g=n_elements(M_g)
  n_s=n_elements(M_s)
  n_tot=n_g+n_s
  if do_bh then begin
    n_bh=n_elements(M_bh)
    n_tot=n_tot+n_bh
  endif
  
  
  ; units conversions here only Mdot_bh in order to have the possibility to select BH particles to plot
  if (do_bh) then begin
    Mdot_bh=Mdot_bh*1e10/hubble*0.98/1e9  
  endif

  if 1 then begin
    device, true=24, retain=2, decomposed=0 
     loadct, 12
     window, win, xsize=700, ysize=700
     plot, [0],[0], XRANGE= [-300,300], YRANGE= [-300,300], psym=3, xstyle=1, ystyle=1, /nodata
     oplot, pos_s[0,*]/(1+z)/hubble, pos_s[1,*]/(1+z)/hubble, psym=3, color=100
     oplot, pos_g[0,*]/(1+z)/hubble, pos_g[1,*]/(1+z)/hubble, psym=3, color=200
    ;if do_bh then oplot, pos_BH[0,where(mdot_bh gt 0.01)], pos_BH[1,where(mdot_bh gt 0.01)], psym=4, color=30
     if do_bh then oplot, pos_BH[0,*], pos_BH[1,*], psym=4, color=cgColor('orange')
     window, win+1, xsize=700, ysize=700
     plot, [0], [0], XRANGE= [-300,300], YRANGE= [-300,300], psym=3, xstyle=1, ystyle=1, /nodata
     oplot, pos_s[0,*]/(1+z)/hubble, pos_s[2,*]/(1+z)/hubble, psym=3, color=100
     oplot, pos_g[0,*]/(1+z)/hubble, pos_g[2,*]/(1+z)/hubble, psym=3, color=200
     ;if do_bh then oplot, pos_BH[0,where(mdot_bh gt 0.01)], pos_BH[2,where(mdot_bh gt 0.01)], psym=4, color=30
     if do_bh then oplot, pos_BH[0,*], pos_BH[2,*], psym=4, color=cgColor('orange')
   endif
 
  ; all units conversions
  pos_g=pos_g/(1+z)/hubble
  pos_s=pos_s/(1+z)/hubble
  if do_bh then begin
    pos_bh=pos_bh/(1+z)/hubble
  endif
  
  n_ele=n_elements(z_g)/n_elements(M_g) ; number of chemical elements followed in the simu
  
  met_g=total(Z_g[1:n_ele-1,*],1)/M_g
  met_s=total(Z_s[1:n_ele-1,*],1)/M_s
  M_g=M_g*1e10/hubble;*cldx_g ; se si usa ind_dusty pensando che solo la massa cold tenga polvere....boh FORSE USARE PER DARE DIRETTAMENTE FRZIONE MOL COME FANNO I MUPPICI?
  M_s=M_s*1e10/hubble
  Mi_s=Mi_s*1e10/hubble
  age_s=galage(z,100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)-galage(1/age_s-1, 100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)
  age_s=age_s*1e-9
  min_allowed_age=1e-4
  jovenes = where(age_s lt min_allowed_age,njov)
  if njov gt 0 then age_s[jovenes]=min_allowed_age
  rho_g=rho_g*(1+z)^3*hubble^2*1e10*1e9;*cldx_g ; se si usa ind_dusty pensando che solo la massa cold tenga polvere....boh
  if do_bh then begin
    M_bh=M_bh*1e10/hubble
  endif
  
 print, keyword_set(write) 
  
  
  if keyword_set(write) then begin ;;;;;Euge
    
    print, 'writing file ', file_out_name
  
    get_lun, uni
    openw, uni, file_out_name
    ;printf, uni, -666.666 , n_tot,  1          ;Euge: imprimo todas las particulas
    printf, uni, -666.666 , n_tot-n_g,  1     ;Euge: No imprimo particulas de gas-> engaño para la desconsideracion  de polvo
    for i=0L, n_s-1 do begin
      printf, uni, format='(i10,9g15.5)', -1, pos_s[*,i], age_s[i], met_s[i], 0.0, Mi_s[i], M_s[i], -666.666
    endfor
   ; Euge: Las siguientes tres lineas debe comentarse cuando no consideramos polvo 
   ; for i=0L, n_g-1 do begin     ; Euge
   ;  printf, uni, format='(i10,9g15.5)', 1, pos_g[*,i], -666.666, met_g[i], rho_g[i], M_g[i], M_g[i], cldx_g[i] ;Euge
   ; endfor ;Euge
    
    if do_bh then begin
      for i=0L, n_bh-1 do begin
        printf, uni, format='(i10,9g15.5)', -1, pos_bh[*,i], 666.666, -666.666, Mdot_bh[i], M_bh[i], M_bh[i], -666.666
      endfor
    endif
    
    close, uni
    free_lun, uni
  
  endif else begin;;;;;Euge

  
  ;  print, 'not writing file';;;;;Euge

  
  endelse;;;;;Euge


;stop
end 



pro many_reg

  for ii = 0, 28 do begin
;  for ii = 25, 27 do begin ; non so perche' cin aveva commentato la 25 in dia_to_g3d
    dia_to_g3d, 0, ii, auto=1, /write
  endfor  
  

end




pro sfr, win, ireg, auto=auto

  ;common COSM, hubble,omega_lambda,omega_m

  ;MIND THE FLAVOR, largheza box!
  
  do_bh=0
  
  box= 2000. ; larghezza fisica del box
  
  simu_dir='/scratch/borgani/Dianoga/'
  regions=['g0016649_G','g0052436_Me14_G','g0144846_Me14_G','g0163178_Me14_G','g0168361_Me14_G','g0272097_G',                $ ; 5
           'g1212639_G','g1483463_G'     ,'g1574117_Me14_G','g1657050_G'     ,'g1680241_G'     ,'g1987669_G',                $ ; 11
           'g2980844_G',                                                                                                     $ ; 12
           'g3327821_G','g3346905_G'     ,'g3888703_G'     ,'g4425770_G'     ,'g4606589_G'     ,'g4915399_G','g5265133_G',   $ ; 19
           'g5503149_G','g5699754_G'     ,'g6287794_G'     ,'g6348555_G'     ,'g6802296_G'     ,                             $ ; 24
           'g7263961_G',                                                                                                     $ ; 25 perche' era commentata??
           'g7358274_G', 'g7570066_G','g7577931_G']   
  flav='csf_bhnew-w500'
  ;flav='csf_w500'
  snap='032'
  
  ;ireg= 2 ; non c'e' snap_032 para el flav bhnew de la regione 1, 3, 4 


  snap_name=simu_dir+regions[ireg]+'/'+flav+'/snap_'+snap
  print, 'SNAP Name: ', snap_name
  sub_name=simu_dir+regions[ireg]+'/'+flav+'/groups_'+snap+'/sub_'+snap
  print, 'SUB  Name: ', sub_name
  
  if keyword_set(auto)  then begin
    auto=auto-1
    readnew, sub_name, spos, 'SPOS'  ; posizioni di tutti i sottoaloni
    readnew, sub_name, fsub, 'FSUB'  ; indici del primo sottoalone di ogni alone negli array di sottoaloni (tipo spos)
    readnew, sub_name, mvir, 'MVIR'  ; masse di ogni alone
    readnew, sub_name, rvir, 'RVIR'  ; masse di ogni alone
    readnew, sub_name, ncon, 'NCON'  ; contaminacion de cada halo, si es distinto de cero esta contaminado 
    
  ;  for i=0,100 do begin
  ;    masse=msub[fsub[i]:fsub[i+1]-1]
  ;    maxmas=max(masse,imax)
  ;    print,  imax
  ;  endfor
  
    mvir = mvir[where(ncon eq 0)]
    fsub = fsub[where(ncon eq 0)]
    
    ind=reverse(sort(mvir))
    fsub=fsub[ind] ; ora fsub e' ordinato in senso decrescente della massa dell'alone
  ;  mvir=mvir[ind] ; non so se servira'
    pos_cdm=spos(*,fsub[auto])
    print, 'pos_cdm', pos_cdm
    
    auto_str=strtrim(string(auto),2)
  endif else begin
    print, 'WARNING !!!! cmd scritto a mano'  
;    pos_cdm=[492777.,498712.,497772.] ; g0016649_G/csf_bhnew-w500 snap 032 altra pos interessante: pos_cdm=[492783.,498704., 497753.]
;    pos_cdm=[497658.,499553.,496432.] ; g0016649_G/csf_w500 snap 041 a ocio perche manca info gropus
;    pos_cdm=[492637.,498740.,497819.] ; g0016649_G/csf_w500 snap 032 a ocio perche manca info gropus molto incerto a z=2
    pos_cdm=[492757.,498711.,497772.] ; g0016649_G/csf_w500 snap 032 passatomi da ste z=2

    auto_str=''
  endelse
  
  box_str=strtrim(string(nint(box)),2)
  if do_bh then begin
    file_out_name='/scratch3/granato/g3d/Input/'+regions[ireg]+'_'+flav+'_'+snap+'_clu'+auto_str+'_L'+box_str+'kpc.dat'    
  endif else begin
    file_out_name='/scratch3/granato/g3d/Input/'+regions[ireg]+'_'+flav+'_'+snap+'_clu'+auto_str+'_L'+box_str+'kpc_noagn.dat'
  endelse  
  
  ;stop
  
  ; test low res cioe' la ris di sempre...
  ;file_out_name="/scratch3/granato/Image/Input/dianoga_test_lr.dat"
  ;snap_name='/scratch/borgani/Dianoga/g0052436_Me14_G/csf_w500/snap_040'
  ;pos_cdm=[499631.,493374.,501427.]
  ;box_half_size=250.
  
  
  ; test low res cioe' la ris di sempre...
  ;file_out_name="/scratch3/granato/g3d/Input/g0016649.dat"
  ;snap_name='/scratch/murante/DianogaHR/g0052436_Me14_G/snap_027'
  ;pos_cdm=[499688.,493465.,501330.] ; file dianoga_test_hr.dat snap 41
  ;pos_cdm=[499971.,498070.,501836.] ; dia2_test_hr.dat snap 038
  ;pos_cdm=[499957.,498125., 497606.] ; g0016649_G/csf_bhnew-w500 snap 032 altra pos interessante: pos_cdm=[492783.,498704., 497753.]
  
  ;box_half_size=125. file dianoga_test_hr.dat snap 41 dia2_test_hr.dat snap 038
  
  
  readnew, snap_name, head, "HEAD"
  z=head.redshift
  hubble=head.hubbleparam
  omega_m=head.omega0
  omega_lambda=head.omegalambda
  
  ; gas
  readnew, snap_name, cldx_g, "CLDX", parttype=0
  readnew, snap_name, pos_g, "POS ", parttype=0
  readnew, snap_name, M_g, "MASS", parttype=0
  readnew, snap_name, Z_g, "Zs  ", parttype=0
  readnew, snap_name, rho_g, "RHO ", parttype=0 ; this is rho_code:  
  readnew, snap_name, T_g, "TEMP", parttype=0
  
  ; rho_phys = rho_code * (1 + z)^3 * m_unit / l_unit^3 * hpar^2
  ; where m_unit=10^10 Msun and l_unit=kpc
  
  ;stars 
  readnew, snap_name, pos_s, "POS ", parttype=4
  readnew, snap_name, M_s, "MASS", parttype=4
  readnew, snap_name, Mi_s, "iM  ", parttype=4
  readnew, snap_name, Z_s, "Zs  ", parttype=4
  readnew, snap_name, age_s, "AGE ", parttype=4
  
  ; BHs
  if do_bh then begin
    readnew, snap_name, pos_bh, "POS ", parttype=5
    readnew, snap_name, M_bh, "BHMA", parttype=5
    ; in snapshot accretion rate given in units of 10^10 h^-1 Msun / (0.98 Gyr). 
    ; 0.98 is the ratio between 3.08 (coming from parsec in cm) and 3.155 (coming from number of seconds in one year)
    readnew, snap_name, Mdot_bh, "BHMD ", parttype=5
    ;readnew, snap_name, sl, "ACRBH", parttype=5
  endif
  
  box_half_gadget=box*(1+z)*hubble/2.
  
  ; select only particles in the desired box
  ind_box_g=where( abs(pos_g[0,*]-pos_cdm[0]) lt box_half_gadget and abs(pos_g[1,*]-pos_cdm[1]) lt box_half_gadget and abs(pos_g[2,*]-pos_cdm[2]) lt box_half_gadget)
  ind_box_s=where( abs(pos_s[0,*]-pos_cdm[0]) lt box_half_gadget and abs(pos_s[1,*]-pos_cdm[1]) lt box_half_gadget and abs(pos_s[2,*]-pos_cdm[2]) lt box_half_gadget)
  pos_g=pos_g[*,ind_box_g] & M_g=M_g[ind_box_g] & Z_g=Z_g[*,ind_box_g] & rho_g=rho_g[ind_box_g] & cldx_g=cldx_g[ind_box_g] & t_g=t_g[ind_box_g]
  pos_s=pos_s[*,ind_box_s] & M_s=M_s[ind_box_s] & Mi_s=Mi_s[ind_box_s] & Z_s=Z_s[*,ind_box_s] & age_s=age_s[ind_box_s] 
  if do_bh then begin
    ind_box_bh=where( abs(pos_bh[0,*]-pos_cdm[0]) lt box_half_gadget and abs(pos_bh[1,*]-pos_cdm[1]) lt box_half_gadget and abs(pos_bh[2,*]-pos_cdm[2]) lt box_half_gadget)
    pos_bh=pos_bh[*,ind_box_bh] & M_bh=M_bh[ind_box_bh] & Mdot_bh=Mdot_bh[ind_box_bh] 
  endif
  
  
  ; select only gas particles which can host dust. cldx_g is the fraction of mass of gas particle in cold phase
  ; if comment this to use all particles, remember also below in units conversion comment multiplications bu cldx_g
  ;ind_dusty=where(cldx_g gt 0 ) ; PRIMA VERSIONE
  ind_dusty=where(cldx_g gt 0 or t_g lt 1e5) ; ALLA MUPPICI..
  pos_g=pos_g[*,ind_dusty] & M_g=M_g[ind_dusty] & Z_g=Z_g[*,ind_dusty] & rho_g=rho_g[ind_dusty] & cldx_g=cldx_g[ind_dusty]
  
  for i=0,2 do begin ; refer coordinates to centre of selected object
    pos_g[i,*]=pos_g[i,*]-pos_cdm[i]
    pos_s[i,*]=pos_s[i,*]-pos_cdm[i]
    if do_bh then begin
      pos_bh[i,*]=pos_bh[i,*]-pos_cdm[i]
    endif
  endfor
  
  ; 
  n_g=n_elements(M_g)
  n_s=n_elements(M_s)
  n_tot=n_g+n_s
  if do_bh then begin
    n_bh=n_elements(M_bh)
    n_tot=n_tot+n_bh
  endif
  
  
    ; all units conversions
;  if (do_bh) then begin
;    Mdot_bh=Mdot_bh*1e10/hubble*0.98/1e9  
;  endif
;  pos_g=pos_g/(1+z)/hubble
;  pos_s=pos_s/(1+z)/hubble
;  if do_bh then begin
;    pos_bh=pos_bh/(1+z)/hubble
;  endif
  
  n_ele=n_elements(z_g)/n_elements(M_g) ; number of chemical elements followed in the simu
  
  met_g=total(Z_g[1:n_ele-1,*],1)/M_g
  met_s=total(Z_s[1:n_ele-1,*],1)/M_s
;  M_g=M_g*1e10/hubble;*cldx_g ; se si usa ind_dusty pensando che solo la massa cold tenga polvere....boh FORSE USARE PER DARE DIRETTAMENTE FRZIONE MOL COME FANNO I MUPPICI?
;  M_s=M_s*1e10/hubble
;  Mi_s=Mi_s*1e10/hubble
;  age_s=galage(z,100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)-galage(1/age_s-1, 100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)
;  age_s=age_s*1e-9
;  min_allowed_age=1e-4
;  age_s[where(age_s lt min_allowed_age)]=min_allowed_age
;  rho_g=rho_g*(1+z)^3*hubble^2*1e10*1e9;*cldx_g ; se si usa ind_dusty pensando che solo la massa cold tenga polvere....boh
;  if do_bh then begin
;    M_bh=M_bh*1e10/hubble
;  endif
 
  dpart=sqrt(pos_s[0,*]^2+pos_s[1,*]^2+pos_s[2,*]^2) 
 
  calc_sfr, age_s, mi_s, 1, dpart, box_half_gadget, 300, z, 5, xbins_s, sfr_s, nxbin_s
  

stop
end 


