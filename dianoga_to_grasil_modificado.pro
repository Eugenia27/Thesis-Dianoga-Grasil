pro dia_to_g3d, win, ireg, biased_progenitor, direct_progenitor, masks, write=write

  n_s    = 0L
  n_g    = 0L
  n_bh   = 0L
  box    = 500.     
  
  simu_dir   = '/gss/gss_work/DRES_murante/CLUSTERS/Dianoga/'
  output_dir = '/pico/scratch/userexternal/cragonef/mferraro/g3d/Input/' 
  print, 'Simu DIR: ', simu_dir

  if ireg EQ 0 then begin
    print, 'D0 DOES NOT EXIST'
    exit
  endif
 
  regions=['none','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10',      $
           'D11','D12','D13','D14','D15','D16','D17','D18','D19','D20',    $
           'D21','D22','D23','D24','D25','D26','D27','D28','D29']

  flav='BH2015'
  snap='032'
  
  snap_name = simu_dir+regions[ireg]+'/'+flav+'/snap_'+snap
  sub_name  = simu_dir+regions[ireg]+'/'+flav+'/Subfind/groups_'+snap+'/sub_'+snap
  print, 'SNAP Name: ', snap_name
  print, 'SUB  Name: ', sub_name

  ; Cosmological SNAP data
  readnew, snap_name, head, "HEAD"
  z           = head.redshift
  hubble      = head.hubbleparam
  omega_m     = head.omega0
  omega_lambda= head.omegalambda
  print, 'REDSHIFT: ', z
 
  ; gas
  readnew, snap_name, cldx_g , "CLDX" , parttype=0
  readnew, snap_name, pos_g  , "POS " , parttype=0
  readnew, snap_name, M_g    , "MASS" , parttype=0
  readnew, snap_name, Z_g    , "Zs  " , parttype=0
  readnew, snap_name, rho_g  , "RHO " , parttype=0  
  readnew, snap_name, T_g    , "TEMP" , parttype=0
 
  ; stars 
  readnew, snap_name, pos_s  , "POS " , parttype=4
  readnew, snap_name, M_s    , "MASS" , parttype=4
  readnew, snap_name, Mi_s   , "iM  " , parttype=4
  readnew, snap_name, Z_s    , "Zs  " , parttype=4
  readnew, snap_name, age_s  , "AGE " , parttype=4
  
  ; BHs
  readnew, snap_name, pos_bh , "POS " , parttype=5
  readnew, snap_name, M_bh   , "BHMA" , parttype=5
  readnew, snap_name, Mdot_bh, "BHMD ", parttype=5
 
  ; SubFind bolcks 
  readnew, sub_name, spos, 'SPOS'  ; posizioni di tutti i sottoaloni
  readnew, sub_name, gpos, 'GPOS'  ; posizioni di tutti i sottoaloni
  readnew, sub_name, smst, 'SMST'  ; masas de cada subhalo
  readnew, sub_name, fsub, 'FSUB'  ; indici del primo sottoalone di ogni alone negli array di sottoaloni (tipo spos)
  readnew, sub_name, ncon, 'NCON'  ; contaminacion de cada halo, si es distinto de cero esta contaminado 
  readnew, sub_name, grnr, 'GRNR'  ; grupo al que pertenece un dado subhalo
  readnew, sub_name, mtot, 'MTOT'  ; masa total de todas las particulas pertenecientes al subhalo
    
  pos_cdm = dblarr(3,2)  
  if biased_progenitor eq 1 then begin
    print,''
    print,'BIASED-PROGENITOR ACTIVATED'
    pos_cdm[*,0] = spos(*,fsub[0])
  endif

  if direct_progenitor eq 1 then begin
    print,''
    print,'DIRECT-PROGENITOR ACTIVATED'
    print,''
    progenitors   = looking_progenitors(ireg, 0, 91, UINT(snap), 1, 1, 0, 0, 1, 0, 15, yhist, mass)
    print,''
    progenitor    = progenitors[1]
    group         = grnr[progenitor]
    main_of_group = fsub[group]

    fname='progenitors.dat'
    openw,21,fname,/append
    printf,21, 'D'+strtrim(string(ireg),2),progenitor,group,main_of_group, FORMAT='(A6,2X,I10,2X,I10,2X,I10)'
    close,21

    if group eq 0 and progenitor eq main_of_group then begin
      print,'DIRECT PROGENITOR IS THE SAME AS THE BIASED PROGENITOR'
    endif else begin
      pos_cdm[*,1] = spos(*,progenitor)
      print,'DIRECT PROGENITOR IS NOT BIASED PROGENITOR'
    endelse 
  endif
    
  for prog=0,1 do begin
    if prog eq 0 and biased_progenitor eq 0 then continue 
    if prog eq 1 and direct_progenitor eq 0 then continue 

    print,'POS CENTER OF MASS OF PROGENITOR', prog,': ', pos_cdm[*,prog]
    progen_type = ''
    if prog eq 0 then progen_type = 'BP'
    if prog eq 1 then progen_type = 'DP'

    box_half_gadget = box*(1+z)*hubble/2.
    file_out_name   = output_dir+regions[ireg]+'_'+snap+'_'+progen_type+'.dat'    

    ; select only particles in the desired box
    if masks eq 0 then begin
      ind_box_g  = where(abs(pos_g[0,*]-pos_cdm[0,prog])  lt box_half_gadget and abs(pos_g[1,*]-pos_cdm[1,prog])  lt box_half_gadget and abs(pos_g[2,*]-pos_cdm[2,prog])  lt box_half_gadget)
      ind_box_s  = where(abs(pos_s[0,*]-pos_cdm[0,prog])  lt box_half_gadget and abs(pos_s[1,*]-pos_cdm[1,prog])  lt box_half_gadget and abs(pos_s[2,*]-pos_cdm[2,prog])  lt box_half_gadget)
      ind_box_bh = where(abs(pos_bh[0,*]-pos_cdm[0,prog]) lt box_half_gadget and abs(pos_bh[1,*]-pos_cdm[1,prog]) lt box_half_gadget and abs(pos_bh[2,*]-pos_cdm[2,prog]) lt box_half_gadget)
      pos_g_  = pos_g[*,ind_box_g]   & M_g_=M_g[ind_box_g]    & Z_g_=Z_g[*,ind_box_g] & rho_g_=rho_g[ind_box_g] & cldx_g_=cldx_g[ind_box_g] & t_g_=t_g[ind_box_g]
      pos_s_  = pos_s[*,ind_box_s]   & M_s_=M_s[ind_box_s]    & Mi_s_=Mi_s[ind_box_s] & Z_s_=Z_s[*,ind_box_s]   & age_s_=age_s[ind_box_s] 
      pos_bh_ = pos_bh[*,ind_box_bh] & M_bh_=M_bh[ind_box_bh] & Mdot_bh_=Mdot_bh[ind_box_bh] 
    endif else begin
    endelse
  
    ind_dusty=where(cldx_g_ gt 0 or t_g_ lt 1e5, nind_dusty)
    print, 'N DUSTY= ', nind_dusty 
    pos_g_ = pos_g_[*,ind_dusty] & M_g_=M_g_[ind_dusty] & Z_g_=Z_g_[*,ind_dusty] & rho_g_=rho_g_[ind_dusty] & cldx_g_=cldx_g[ind_dusty]
  
    for i=0,2 do begin                   ; refer coordinates to centre of selected object
      pos_g_[i,*]  = pos_g_[i,*] -pos_cdm[i,prog]
      pos_s_[i,*]  = pos_s_[i,*] -pos_cdm[i,prog]
      pos_bh_[i,*] = pos_bh_[i,*]-pos_cdm[i,prog]
    endfor
   
    n_g   = n_elements(M_g_)
    n_s   = n_elements(M_s_)
    n_bh  = n_elements(M_bh_)
    n_tot = n_g+n_s+n_bh
  
    ; units conversions here only Mdot_bh in order to have the possibility to select BH particles to plot
    Mdot_bh_ = Mdot_bh_*1e10/hubble*0.98/1e9  

    if 1 then begin
      device, true=24, retain=2, decomposed=0 
       loadct, 12
       window, win, xsize=700, ysize=700
       plot,[0],[0], XRANGE= [-300,300], YRANGE= [-300,300], psym=3, xstyle=1, ystyle=1, /nodata
       oplot,pos_s_[0,*]/(1+z)/hubble, pos_s_[1,*]/(1+z)/hubble, psym=3, color=100
       oplot,pos_g_[0,*]/(1+z)/hubble, pos_g_[1,*]/(1+z)/hubble, psym=3, color=200
       oplot,pos_BH_[0,*], pos_BH[1,*], psym=4, color=cgColor('orange')
       window, win+1, xsize=700, ysize=700
       plot, [0], [0], XRANGE= [-300,300], YRANGE= [-300,300], psym=3, xstyle=1, ystyle=1, /nodata
       oplot, pos_s_[0,*]/(1+z)/hubble, pos_s_[2,*]/(1+z)/hubble, psym=3, color=100
       oplot, pos_g_[0,*]/(1+z)/hubble, pos_g_[2,*]/(1+z)/hubble, psym=3, color=200
       oplot, pos_BH_[0,*], pos_BH[2,*], psym=4, color=cgColor('orange')
     endif
 
    ; all units conversions
    pos_g_  = pos_g_/(1+z)/hubble
    pos_s_  = pos_s_/(1+z)/hubble
    pos_bh_ = pos_bh_/(1+z)/hubble
    n_ele = n_elements(z_g_)/n_elements(M_g_)        ; number of chemical elements followed in the simulation
  
    met_g = total(Z_g_[1:n_ele-1,*],1)/M_g_
    met_s = total(Z_s_[1:n_ele-1,*],1)/M_s_
    M_g_   = M_g_*1e10/hubble                          
    M_s_   = M_s_*1e10/hubble
    Mi_s_  = Mi_s_*1e10/hubble
    age_s_ = galage(z,100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)-galage(1/age_s_-1, 100,H0=hubble*100,lambda0=omega_lambda,omega_m=omega_m)
    age_s_ = age_s*1e-9
    min_allowed_age = 1e-4
    jovenes = where(age_s_ lt min_allowed_age,njov)
    if njov gt 0 then age_s_[jovenes] = min_allowed_age
    rho_g_ = rho_g_*(1+z)^3*hubble^2*1e10*1e9        
    M_bh_  = M_bh_*1e10/hubble
  
    if keyword_set(write) then begin
      print, 'WRITING FILE ', file_out_name
      get_lun, uni
      openw, uni, file_out_name
      ;printf, uni, -666.666 , n_tot,  1        ;Guardo todas las particulas
      printf, uni, -666.666 , n_tot-n_g,  1     ;No guardo particulas de gas-> engaño para desconsiderar polvo
      for i=0L, n_s-1 do begin
        printf, uni, format='(i10,9g15.5)', -1, pos_s_[*,i], age_s_[i], met_s[i], 0.0, Mi_s_[i], M_s_[i], -666.666
      endfor
      ; No polvo: comentar siguiente blucle for 
      ; for i=0L, n_g-1 do begin     
      ;  printf, uni, format='(i10,9g15.5)', 1, pos_g_[*,i], -666.666, met_g[i], rho_g_[i], M_g_[i], M_g_[i], cldx_g_[i] 
      ; endfor 
      for i=0L, n_bh-1 do begin
        printf, uni, format='(i10,9g15.5)', -1, pos_bh_[*,i], 666.666, -666.666, Mdot_bh_[i], M_bh_[i], M_bh_[i], -666.666
      endfor
      close, uni
      free_lun, uni
    endif
    n_g   = 0L 
    n_s   = 0L
    n_bh  = 0L
    n_tot = 0L
    nind_dusty    = 0
    ind_box_g     = 0
    ind_box_s     = 0
    ind_box_bh    = 0
    ind_box_dusty = 0
    pos_g_ = [ ]
    pos_s_ = [ ]
    pos_bh_= [ ]
    M_g_   = [ ]                         
    M_s_   = [ ]
    M_bh_  = [ ]
    Z_g_   = [ ]
    T_g_   = [ ]
    Z_s_   = [ ]
    rho_g_ = [ ]
    Mi_s_  = [ ]
    age_s_ = [ ]
    jovenes= [ ]
    Modot_bh_= [ ]
    njov   = 0 
    n_ele  = 0
    met_g  = 0.
    met_s  = 0.
  endfor
end


